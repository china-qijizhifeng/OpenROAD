// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Priority queue of the vertex search (PathSearch). No drt dependencies.
//
// A 4-ary min-heap whose groups of four children start on a cache line
// (the root sits at index 3 of a 64-byte aligned array, so the children of
// a node are the indices 4 (k + 1) .. 4 (k + 1) + 3), and the last pushed
// item held outside the heap: it enters the heap only when another item is
// pushed before it is popped. pop() returns the minimum of the heap and the
// held item, so the pop sequence is that of a plain heap for any strict
// total order Item::operator>; a run label whose continuation is the next
// pop costs no heap operation.

#pragma once

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <new>
#include <vector>

namespace drt::br {

template <typename T>
struct CacheAlignedAllocator
{
  using value_type = T;
  static constexpr std::align_val_t kAlignment{64};
  CacheAlignedAllocator() = default;
  template <typename U>
  CacheAlignedAllocator(const CacheAlignedAllocator<U>&)  // NOLINT
  {
  }
  T* allocate(const size_t n)
  {
    return static_cast<T*>(::operator new(n * sizeof(T), kAlignment));
  }
  void deallocate(T* p, size_t) { ::operator delete(p, kAlignment); }
  template <typename U>
  bool operator==(const CacheAlignedAllocator<U>&) const
  {
    return true;
  }
};

// (z, i, j, dir) of a vertex-search item in one word whose order is the
// lexicographic order of (z, i, j, dir): z above i above j above two bits
// of dir + 1 (dir in -1 .. 2). Field widths follow the largest layer,
// track and position counts.
class VertexOrder
{
 public:
  // False when the widths exceed 64 bits.
  bool setup(const int layers, const int tracks, const int positions)
  {
    const int bits_j = std::bit_width(static_cast<unsigned>(std::max(positions, 1)));
    const int bits_i = std::bit_width(static_cast<unsigned>(std::max(tracks, 1)));
    const int bits_z = std::bit_width(static_cast<unsigned>(std::max(layers, 1)));
    mask_j_ = (uint64_t{1} << bits_j) - 1;
    mask_i_ = (uint64_t{1} << bits_i) - 1;
    shift_i_ = 2 + bits_j;
    shift_z_ = shift_i_ + bits_i;
    return shift_z_ + bits_z <= 64;
  }
  uint64_t pack(const int z, const int i, const int j, const int dir) const
  {
    return (static_cast<uint64_t>(z) << shift_z_)
           | (static_cast<uint64_t>(i) << shift_i_)
           | (static_cast<uint64_t>(j) << 2) | static_cast<uint64_t>(dir + 1);
  }
  int z(const uint64_t order) const
  {
    return static_cast<int>(order >> shift_z_);
  }
  int i(const uint64_t order) const
  {
    return static_cast<int>((order >> shift_i_) & mask_i_);
  }
  int j(const uint64_t order) const
  {
    return static_cast<int>((order >> 2) & mask_j_);
  }
  static int dir(const uint64_t order)
  {
    return static_cast<int>(order & 3) - 1;
  }

 private:
  int shift_i_ = 2;
  int shift_z_ = 2;
  uint64_t mask_i_ = 0;
  uint64_t mask_j_ = 0;
};

// Item: trivially copyable with a strict total order operator> (no two
// queued items compare equal).
template <typename Item>
class SearchQueue
{
 public:
  void clear()
  {
    heap_.resize(kRoot);
    held_ = false;
  }
  bool empty() const { return !held_ && heap_.size() == kRoot; }
  size_t size() const { return heap_.size() - kRoot + (held_ ? 1 : 0); }
  void push(const Item& item)
  {
    if (held_) {
      heapPush(last_);
    }
    last_ = item;
    held_ = true;
  }
  // The minimum; the queue must not be empty.
  Item pop()
  {
    if (held_ && (heap_.size() == kRoot || heap_[kRoot] > last_)) {
      held_ = false;
      held_pops_++;
      return last_;
    }
    return heapPop();
  }
  // Pops answered by the held item (no heap operation), since construction.
  int64_t heldPops() const { return held_pops_; }

 private:
  static constexpr size_t kRoot = 3;

  void heapPush(const Item& item)
  {
    heap_.push_back(item);
    Item* h = heap_.data() + kRoot;
    size_t k = heap_.size() - 1 - kRoot;
    while (k > 0) {
      const size_t parent = (k - 1) / 4;
      if (!(h[parent] > item)) {
        break;
      }
      h[k] = h[parent];
      k = parent;
    }
    h[k] = item;
  }
  Item heapPop()
  {
    Item* h = heap_.data() + kRoot;
    const Item top = h[0];
    const Item last = heap_.back();
    heap_.pop_back();
    const size_t n = heap_.size() - kRoot;
    if (n == 0) {
      return top;
    }
    size_t k = 0;
    while (true) {
      const size_t first = 4 * k + 1;
      if (first >= n) {
        break;
      }
      size_t best = first;
      const size_t end = first + 4 < n ? first + 4 : n;
      for (size_t c = first + 1; c < end; c++) {
        if (h[best] > h[c]) {
          best = c;
        }
      }
      if (!(last > h[best])) {
        break;
      }
      h[k] = h[best];
      k = best;
    }
    h[k] = last;
    return top;
  }

  std::vector<Item, CacheAlignedAllocator<Item>> heap_
      = std::vector<Item, CacheAlignedAllocator<Item>>(kRoot);
  Item last_{};
  bool held_ = false;
  int64_t held_pops_ = 0;
};

}  // namespace drt::br
