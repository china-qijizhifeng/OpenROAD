// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2018-2025, The OpenROAD Authors

#pragma once

#include <cstddef>
#include <utility>
#include <vector>

namespace gpl {

class FloatPoint
{
 public:
  float x = 0;
  float y = 0;
  FloatPoint() = default;
  FloatPoint(float x, float y) : x(x), y(y) {}
};

// Optimizer coordinates and gradients own separate x/y columns throughout
// backtracking and rollback. Scalar row views are only an API convenience;
// kernels can access x()/y() directly without repacking.
class PointColumns
{
 public:
  struct Reference
  {
    float& x;
    float& y;
    operator FloatPoint() const { return {x, y}; }
    Reference& operator=(FloatPoint point)
    {
      x = point.x;
      y = point.y;
      return *this;
    }
    Reference& operator=(const Reference& other)
    {
      return *this = static_cast<FloatPoint>(other);
    }
  };

  size_t size() const { return x_.size(); }
  bool empty() const { return x_.empty(); }
  void resize(size_t size, FloatPoint value = {})
  {
    x_.resize(size, value.x);
    y_.resize(size, value.y);
  }
  void emplace_back(FloatPoint value = {})
  {
    x_.push_back(value.x);
    y_.push_back(value.y);
  }
  void pop_back()
  {
    x_.pop_back();
    y_.pop_back();
  }
  void swapRows(size_t a, size_t b)
  {
    std::swap(x_[a], x_[b]);
    std::swap(y_[a], y_[b]);
  }
  Reference operator[](size_t i) { return {x_[i], y_[i]}; }
  FloatPoint operator[](size_t i) const { return {x_[i], y_[i]}; }
  std::vector<float>& x() { return x_; }
  std::vector<float>& y() { return y_; }
  const std::vector<float>& x() const { return x_; }
  const std::vector<float>& y() const { return y_; }

 private:
  std::vector<float> x_, y_;
};

}  // namespace gpl
