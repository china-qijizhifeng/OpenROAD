// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors
#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

namespace gpl {

// Column storage owns numerical state. Row handles keep a slot, never a pointer
// into a column, so a topology-boundary reallocation cannot invalidate handles.
// Allocation/release is restricted to the serial topology update phase.
// Parallel numerical kernels may access distinct rows or use their existing
// atomics.
template <typename... Fields>
class SoaStorage
{
 public:
  size_t allocate()
  {
    if (!free_.empty()) {
      const size_t slot = free_.back();
      free_.pop_back();
      return slot;
    }
    const size_t slot = std::get<0>(columns_).size();
    // Releases run in destructors: provision their freelist storage here.
    if (free_.capacity() <= slot) {
      free_.reserve(std::max<size_t>(16, (slot + 1) * 2));
    }
    std::apply([](auto&... column) { (column.emplace_back(), ...); }, columns_);
    return slot;
  }

  void release(size_t slot) { free_.push_back(slot); }

  template <size_t Field>
  auto& at(size_t slot)
  {
    return std::get<Field>(columns_)[slot];
  }

  template <size_t Field>
  const auto& at(size_t slot) const
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

  void copy(size_t to, size_t from)
  {
    std::apply([&](auto&... column) { ((column[to] = column[from]), ...); },
               columns_);
  }

  void clearRow(size_t slot)
  {
    std::apply([&](auto&... column) { ((column[slot] = {}), ...); }, columns_);
  }

 private:
  std::tuple<std::vector<Fields>...> columns_;
  std::vector<size_t> free_;
};

// Value semantics matter: filler snapshots and swap/remove operations copy
// objects. A copy gets its own numerical row; a move transfers the row. Shared
// ownership lets a removed filler snapshot outlive the owning region's vector.
template <typename Store>
class SoaRow
{
 public:
  explicit SoaRow(std::shared_ptr<Store> store)
      : store_(std::move(store)), slot_(store_->allocate())
  {
    store_->clearRow(slot_);
  }

  SoaRow(const SoaRow& other) : SoaRow(other.store_)
  {
    store_->copy(slot_, other.slot_);
  }

  SoaRow(SoaRow&& other) noexcept
      : store_(std::move(other.store_)), slot_(other.slot_)
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
    if (store_) {
      store_->release(slot_);
    }
  }

  template <size_t Field>
  auto& get()
  {
    return store_->template at<Field>(slot_);
  }

  template <size_t Field>
  const auto& get() const
  {
    return std::as_const(*store_).template at<Field>(slot_);
  }

  Store& storage() { return *store_; }
  const Store& storage() const { return *store_; }
  size_t slot() const { return slot_; }

 private:
  void swap(SoaRow& other) noexcept
  {
    store_.swap(other.store_);
    std::swap(slot_, other.slot_);
  }

  std::shared_ptr<Store> store_;
  size_t slot_;
};

}  // namespace gpl
