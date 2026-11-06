// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors
#pragma once

#include <algorithm>
#include <cassert>
#include <compare>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iterator>
#include <limits>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace gpl {

struct SoaViewTag
{
};
template <typename Object>
class SoaRef;
template <typename Object>
class SoaRows;
template <typename Object>
class SoaPool;

template <typename Object>
class SoaRef
{
 public:
  SoaRef() = default;
  SoaRef(std::nullptr_t) {}
  SoaRef(void* storage, size_t row) : storage_(storage), row_(row) {}
  struct Arrow
  {
    Object view;
    Object* operator->() { return std::addressof(view); }
  };
  Arrow operator->() const { return {operator*()}; }
  Object operator*() const { return Object(SoaViewTag{}, storage_, row_); }
  explicit operator bool() const { return storage_ != nullptr; }
  void* storage() const { return storage_; }
  size_t row() const { return row_; }
  bool operator==(const SoaRef&) const = default;
  friend size_t hash_value(SoaRef value)
  {
    return std::hash<void*>{}(value.storage_)
           ^ (std::hash<size_t>{}(value.row_) << 1);
  }
  auto operator<=>(const SoaRef&) const = default;

 private:
  void* storage_ = nullptr;
  size_t row_ = 0;
};

// Column-wide table identity: every relationship stores just a 32-bit row.
// Temporary proxies provide the reference API without a table pointer per pin.
template <typename Object>
class SoaReferenceColumn
{
 public:
  using value_type = SoaRef<Object>;
  class Reference : public value_type
  {
   public:
    Reference(SoaReferenceColumn& column, size_t index)
        : value_type(std::as_const(column)[index]),
          column_(&column),
          index_(index)
    {
    }
    Reference& operator=(value_type value)
    {
      column_->set(index_, value);
      static_cast<value_type&>(*this) = value;
      return *this;
    }
    Reference& operator=(const Reference& other)
    {
      return *this = static_cast<value_type>(other);
    }

   private:
    SoaReferenceColumn* column_;
    size_t index_;
  };
  // Bind at the serial topology boundary before distinct-row parallel writes.
  void bindTable(void* storage)
  {
    if (storage_ != nullptr && storage_ != storage) {
      throw std::logic_error("GPL relation spans different SoA tables");
    }
    storage_ = storage;
  }
  size_t size() const { return rows_.size(); }
  size_t capacity() const { return rows_.capacity(); }
  size_t allocatedBytes() const { return rows_.capacity() * sizeof(uint32_t); }
  void reserve(size_t count) { rows_.reserve(count); }
  void emplace_back() { rows_.push_back(kNoRow); }
  Reference operator[](size_t index) { return {*this, index}; }
  value_type operator[](size_t index) const
  {
    return rows_[index] == kNoRow ? value_type{}
                                  : value_type{storage_, rows_[index]};
  }

 private:
  void set(size_t index, value_type value)
  {
    if (!value) {
      rows_[index] = kNoRow;
      return;
    }
    if (storage_ == nullptr) {
      storage_ = value.storage();
    }
    if (storage_ != value.storage()) {
      throw std::logic_error("GPL relation spans different SoA tables");
    }
    if (value.row() >= kNoRow) {
      throw std::length_error("GPL SoA row limit exceeded");
    }
    rows_[index] = value.row();
  }
  static constexpr uint32_t kNoRow = std::numeric_limits<uint32_t>::max();
  void* storage_ = nullptr;
  std::vector<uint32_t> rows_;
};

template <typename Field>
struct SoaColumnType
{
  using type = std::vector<Field>;
};
template <typename Object>
struct SoaColumnType<SoaRef<Object>>
{
  using type = SoaReferenceColumn<Object>;
};

// Column storage owns numerical state. Row handles keep a slot, never a pointer
// into a column, so a topology-boundary reallocation cannot invalidate handles.
// Allocation/release is restricted to the serial topology update phase.
// Parallel numerical kernels may access distinct rows or use their existing
// atomics.
template <typename... Fields>
class SoaStorage : public std::enable_shared_from_this<SoaStorage<Fields...>>
{
 public:
  void reserveRows(size_t count)
  {
    if (count > std::numeric_limits<uint32_t>::max()) {
      throw std::length_error("GPL SoA row limit exceeded");
    }
    std::apply([&](auto&... column) { (column.reserve(count), ...); },
               columns_);
    free_.reserve(count);
  }

  void reserveAdditional(size_t count)
  {
    const size_t extra = count > free_.size() ? count - free_.size() : 0;
    reserveRows(std::get<0>(columns_).size() + extra);
  }

  size_t allocate()
  {
    if (!free_.empty()) {
      const size_t slot = free_.back();
      free_.pop_back();
      return slot;
    }
    const size_t slot = std::get<0>(columns_).size();
    if (slot >= std::numeric_limits<uint32_t>::max()) {
      throw std::length_error("GPL SoA row limit exceeded");
    }
    // Known batches reserve once. Unpredicted callbacks grow all columns
    // together, including destructor-safe freelist capacity.
    if (std::get<0>(columns_).capacity() <= slot) {
      reserveRows(std::min<size_t>(
          std::numeric_limits<uint32_t>::max(),
          std::max<size_t>(16, slot + std::max<size_t>(1, slot / 2))));
    }
    std::apply([](auto&... column) { (column.emplace_back(), ...); }, columns_);
    return slot;
  }

  void release(size_t slot)
  {
    clearRow(slot);
    free_.push_back(static_cast<uint32_t>(slot));
  }

  template <size_t Field>
  decltype(auto) at(size_t slot)
  {
    return std::get<Field>(columns_)[slot];
  }

  template <size_t Field>
  decltype(auto) at(size_t slot) const
  {
    return std::get<Field>(columns_)[slot];
  }

  template <size_t Field>
  auto& column()
  {
    return std::get<Field>(columns_);
  }

  template <size_t Field>
  const auto& column() const
  {
    return std::get<Field>(columns_);
  }

  size_t size() const { return std::get<0>(columns_).size(); }
  size_t liveSize() const { return size() - free_.size(); }
  size_t capacity() const { return std::get<0>(columns_).capacity(); }
  size_t allocatedBytes() const
  {
    size_t bytes = free_.capacity() * sizeof(uint32_t);
    auto count = [&]<typename Column>(const Column& column) {
      if constexpr (requires { column.allocatedBytes(); }) {
        bytes += column.allocatedBytes();
      } else {
        bytes += column.capacity() * sizeof(typename Column::value_type);
        if constexpr (requires { column[0].allocatedBytes(); }) {
          for (const auto& value : column) {
            bytes += value.allocatedBytes();
          }
        }
      }
    };
    std::apply([&](const auto&... column) { (count(column), ...); }, columns_);
    return bytes;
  }

  void copy(size_t to, size_t from)
  {
    std::apply([&](auto&... column) { ((column[to] = column[from]), ...); },
               columns_);
  }

  void clearRow(size_t slot)
  {
    std::apply(
        [&](auto&... column) {
          ((column[slot]
            = typename std::remove_reference_t<decltype(column)>::value_type{}),
           ...);
        },
        columns_);
  }

 private:
  std::tuple<typename SoaColumnType<Fields>::type...> columns_;
  std::vector<uint32_t> free_;
};

// Value semantics matter: filler snapshots and swap/remove operations copy
// objects. A copy gets its own numerical row; a move transfers the row. Shared
// ownership lets a removed filler snapshot outlive the owning region's vector.
template <typename Store>
class SoaRow
{
 public:
  SoaRow(SoaViewTag, void* store, size_t slot)
      : store_(std::shared_ptr<Store>{}, static_cast<Store*>(store)),
        slot_(slot),
        owns_(false)
  {
  }

  explicit SoaRow(std::shared_ptr<Store> store)
      : store_(std::move(store)), slot_(store_->allocate())
  {
    store_->clearRow(slot_);
  }

  SoaRow(const SoaRow& other) : SoaRow(other.owner())
  {
    store_->copy(slot_, other.slot_);
  }

  SoaRow(SoaRow&& other) noexcept
      : store_(std::move(other.store_)), slot_(other.slot_), owns_(other.owns_)
  {
  }

  SoaRow& operator=(const SoaRow& other)
  {
    if (this != &other) {
      SoaRow copy(other);
      swap(copy);
    }
    return *this;
  }

  SoaRow& operator=(SoaRow&& other) noexcept
  {
    if (this != &other) {
      SoaRow moved(std::move(other));
      swap(moved);
    }
    return *this;
  }

  ~SoaRow()
  {
    if (store_ && owns_) {
      store_->release(slot_);
    }
  }

  template <size_t Field>
  decltype(auto) get()
  {
    return store_->template at<Field>(slot_);
  }

  template <size_t Field>
  decltype(auto) get() const
  {
    return std::as_const(*store_).template at<Field>(slot_);
  }

  Store& storage() { return *store_; }
  const Store& storage() const { return *store_; }
  size_t slot() const { return slot_; }
  Store* storageIdentity() const { return store_.get(); }
  std::shared_ptr<Store> owner() const
  {
    return std::static_pointer_cast<Store>(store_->shared_from_this());
  }
  bool owns() const { return owns_; }
  void disown() { owns_ = false; }

 private:
  void swap(SoaRow& other) noexcept
  {
    store_.swap(other.store_);
    std::swap(slot_, other.slot_);
    std::swap(owns_, other.owns_);
  }

  std::shared_ptr<Store> store_;
  size_t slot_;
  bool owns_ = true;
};

// References identify rows, never addresses of objects or column elements.
// The owning placement table keeps storage alive across topology callbacks.

// Dereference returns a borrowed value view. No reference points inside a
// temporary iterator, including OpenMP's indexed range-loop lowering.
template <typename Range, typename Value>
class SoaIterator
{
 public:
  using difference_type = std::ptrdiff_t;
  using value_type = Value;
  using reference = Value;
  using pointer = void;
  using iterator_category = std::random_access_iterator_tag;
  SoaIterator() = default;
  SoaIterator(const Range* range, difference_type index)
      : range_(range), index_(index)
  {
  }
  Value operator*() const { return (*range_)[index_]; }
  Value operator[](difference_type offset) const
  {
    return (*range_)[index_ + offset];
  }
  SoaIterator& operator+=(difference_type n)
  {
    index_ += n;
    return *this;
  }
  SoaIterator& operator-=(difference_type n) { return *this += -n; }
  SoaIterator& operator++() { return *this += 1; }
  SoaIterator& operator--() { return *this -= 1; }
  SoaIterator operator++(int)
  {
    auto old = *this;
    ++*this;
    return old;
  }
  SoaIterator operator--(int)
  {
    auto old = *this;
    --*this;
    return old;
  }
  SoaIterator operator+(difference_type n) const
  {
    auto it = *this;
    return it += n;
  }
  SoaIterator operator-(difference_type n) const
  {
    auto it = *this;
    return it -= n;
  }
  friend SoaIterator operator+(difference_type n, SoaIterator it)
  {
    return it += n;
  }
  difference_type operator-(const SoaIterator& other) const
  {
    return index_ - other.index_;
  }
  bool operator==(const SoaIterator& other) const
  {
    return range_ == other.range_ && index_ == other.index_;
  }
  auto operator<=>(const SoaIterator& other) const
  {
    return index_ <=> other.index_;
  }

 private:
  const Range* range_ = nullptr;
  difference_type index_ = 0;
};

// Homogeneous adjacency stores only row IDs. The table identity is shared by
// the whole list; there are no per-edge pointers or persistent object handles.
template <typename Object>
class SoaRows
{
 public:
  using value_type = SoaRef<Object>;
  using iterator = SoaIterator<SoaRows, value_type>;
  SoaRows() = default;
  SoaRows(std::initializer_list<value_type> values)
  {
    reserve(values.size());
    for (value_type value : values) {
      push_back(value);
    }
  }
  size_t size() const { return rows_.size(); }
  size_t allocatedBytes() const { return rows_.capacity() * sizeof(uint32_t); }
  bool empty() const { return rows_.empty(); }
  void reserve(size_t count) { rows_.reserve(count); }
  void clear() { rows_.clear(); }
  value_type operator[](size_t index) const
  {
    return rows_[index] == std::numeric_limits<uint32_t>::max()
               ? value_type{}
               : value_type{storage_, rows_[index]};
  }
  void resize(size_t count, value_type value = {})
  {
    bind(value);
    rows_.resize(count,
                 value ? static_cast<uint32_t>(value.row())
                       : std::numeric_limits<uint32_t>::max());
  }
  void pop_back() { rows_.pop_back(); }
  void set(size_t index, value_type value)
  {
    if (!value) {
      rows_[index] = std::numeric_limits<uint32_t>::max();
      return;
    }
    bind(value);
    rows_[index] = static_cast<uint32_t>(value.row());
  }
  value_type front() const { return (*this)[0]; }
  value_type back() const { return (*this)[size() - 1]; }
  iterator begin() const { return {this, 0}; }
  iterator end() const { return {this, static_cast<std::ptrdiff_t>(size())}; }
  void push_back(value_type value)
  {
    bind(value);
    rows_.push_back(value ? static_cast<uint32_t>(value.row())
                          : std::numeric_limits<uint32_t>::max());
  }
  void emplace_back(value_type value) { push_back(value); }
  const std::vector<uint32_t>& rows() const { return rows_; }
  void* storage() const { return storage_; }

 private:
  void bind(value_type value)
  {
    if (!value) {
      return;
    }
    if (storage_ == nullptr) {
      storage_ = value.storage();
    }
    if (storage_ != value.storage()) {
      throw std::logic_error("GPL adjacency spans different SoA tables");
    }
    if (value.row() >= std::numeric_limits<uint32_t>::max()) {
      throw std::length_error("GPL SoA row limit exceeded");
    }
  }
  void* storage_ = nullptr;
  std::vector<uint32_t> rows_;
};

// Owns live rows and their iteration order, rather than one C++ object per row.
// Swap/remove permutes IDs only, so all external row references remain stable.
template <typename Object>
class SoaPool
{
 public:
  using iterator = SoaIterator<SoaPool, Object>;
  using Store = typename Object::Store;
  SoaPool() = default;
  explicit SoaPool(std::shared_ptr<Store> store) : store_(std::move(store)) {}
  SoaPool(const SoaPool& other) : store_(other.store_)
  {
    reserve(other.size());
    for (const auto& view : other) {
      push_back(view);
    }
  }
  SoaPool(SoaPool&&) noexcept = default;
  SoaPool& operator=(SoaPool other) noexcept
  {
    swap(other);
    return *this;
  }
  ~SoaPool() { clear(); }
  void swap(SoaPool& other) noexcept
  {
    store_.swap(other.store_);
    rows_.swap(other.rows_);
    positions_.swap(other.positions_);
  }
  size_t size() const { return rows_.size(); }
  size_t allocatedBytes() const
  {
    return (rows_.capacity() + positions_.capacity()) * sizeof(uint32_t);
  }
  bool empty() const { return rows_.empty(); }
  SoaRef<Object> ref(size_t index) const
  {
    return {store_.get(), rows_[index]};
  }
  Object operator[](size_t index) const
  {
    return Object(SoaViewTag{}, store_.get(), rows_[index]);
  }
  Object front() const { return (*this)[0]; }
  Object back() const { return (*this)[size() - 1]; }
  iterator begin() const { return {this, 0}; }
  iterator end() const { return {this, static_cast<std::ptrdiff_t>(size())}; }
  void reserve(size_t count)
  {
    rows_.reserve(count);
    if (!positions_.empty()) {
      positions_.reserve(count);
    }
    if (store_ && count > size()) {
      store_->reserveAdditional(count - size());
    }
  }
  void push_back(const Object& object)
  {
    Object copy(object);
    push_back(std::move(copy));
  }
  void push_back(Object&& object)
  {
    auto& row = object.numericRow();
    if (!row.owns()) {
      Object copy(object);
      push_back(std::move(copy));
      return;
    }
    if (!store_) {
      store_ = row.owner();
    }
    assert(store_.get() == &row.storage());
    appendRow(row.slot());
    row.disown();
  }
  template <typename... Args>
  void emplace_back(Args&&... args)
  {
    Object object(std::forward<Args>(args)...);
    push_back(std::move(object));
  }
  void reorderRows(std::vector<uint32_t> rows)
  {
    assert(rows.size() == rows_.size());
    rows_ = std::move(rows);
    positions_.clear();
  }
  void swapRows(size_t a, size_t b)
  {
    std::swap(rows_[a], rows_[b]);
    if (!positions_.empty()) {
      positions_[rows_[a]] = a;
      positions_[rows_[b]] = b;
    }
  }
  void pop_back()
  {
    if (!positions_.empty()) {
      positions_[rows_.back()] = std::numeric_limits<uint32_t>::max();
    }
    store_->release(rows_.back());
    rows_.pop_back();
  }
  void clear()
  {
    while (!empty()) {
      pop_back();
    }
  }
  size_t indexOf(SoaRef<Object> ref) const
  {
    assert(ref.storage() == store_.get());
    if (positions_.empty() && !rows_.empty()) {
      positions_.assign(
          static_cast<size_t>(*std::max_element(rows_.begin(), rows_.end()))
              + 1,
          std::numeric_limits<uint32_t>::max());
      for (size_t i = 0; i < size(); ++i) {
        positions_[rows_[i]] = i;
      }
    }
    assert(ref.row() < positions_.size());
    return positions_[ref.row()];
  }
  const std::vector<uint32_t>& rows() const { return rows_; }

 private:
  void appendRow(size_t row)
  {
    if (row >= std::numeric_limits<uint32_t>::max()) {
      throw std::length_error("GPL SoA row limit exceeded");
    }
    if (!positions_.empty()) {
      if (positions_.size() <= row) {
        positions_.resize(row + 1, std::numeric_limits<uint32_t>::max());
      }
      positions_[row] = rows_.size();
    }
    rows_.push_back(static_cast<uint32_t>(row));
  }
  std::shared_ptr<Store> store_;
  std::vector<uint32_t> rows_;
  mutable std::vector<uint32_t> positions_;
};

}  // namespace gpl

namespace std {
template <typename T>
struct hash<gpl::SoaRef<T>>
{
  size_t operator()(gpl::SoaRef<T> ref) const
  {
    return hash<void*>{}(ref.storage()) ^ (hash<size_t>{}(ref.row()) << 1);
  }
};
}  // namespace std
