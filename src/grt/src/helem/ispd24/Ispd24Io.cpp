// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "Ispd24Io.h"

#include <cctype>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "core/GCellGraph.h"

namespace grt::helem {

namespace {

bool readWholeFile(const std::string& path,
                   std::string& content,
                   std::string& error)
{
  FILE* file = std::fopen(path.c_str(), "rb");
  if (file == nullptr) {
    error = "cannot open " + path + ": " + std::strerror(errno);
    return false;
  }
  std::fseek(file, 0, SEEK_END);
  const long size = std::ftell(file);
  std::fseek(file, 0, SEEK_SET);
  content.resize(size);
  const size_t got = size > 0 ? std::fread(content.data(), 1, size, file) : 0;
  std::fclose(file);
  if (got != static_cast<size_t>(size)) {
    error = "short read on " + path;
    return false;
  }
  return true;
}

// Minimal cursor over a text buffer.
class Cursor
{
 public:
  explicit Cursor(const std::string& text)
      : position_(text.data()), end_(text.data() + text.size())
  {
  }

  bool atEnd() const { return position_ >= end_; }

  void skipSpace()
  {
    while (position_ < end_
           && std::isspace(static_cast<unsigned char>(*position_))) {
      ++position_;
    }
  }

  bool readDouble(double& value)
  {
    skipSpace();
    if (position_ >= end_) {
      return false;
    }
    char* next = nullptr;
    value = std::strtod(position_, &next);
    if (next == position_) {
      return false;
    }
    position_ = next;
    return true;
  }

  bool readInt64(int64_t& value)
  {
    skipSpace();
    if (position_ >= end_) {
      return false;
    }
    char* next = nullptr;
    value = std::strtoll(position_, &next, 10);
    if (next == position_) {
      return false;
    }
    position_ = next;
    return true;
  }

  // Returns the rest of the current line and moves past its newline.
  std::string readLine()
  {
    const char* start = position_;
    while (position_ < end_ && *position_ != '\n') {
      ++position_;
    }
    std::string line(start, position_);
    if (position_ < end_) {
      ++position_;
    }
    if (!line.empty() && line.back() == '\r') {
      line.pop_back();
    }
    return line;
  }

  void skipLineRemainder()
  {
    while (position_ < end_ && *position_ != '\n') {
      ++position_;
    }
    if (position_ < end_) {
      ++position_;
    }
  }

  const char* position() const { return position_; }
  const char* end() const { return end_; }
  void setPosition(const char* position) { position_ = position; }

 private:
  const char* position_;
  const char* end_;
};

}  // namespace

bool readIspd24Resource(const std::string& path,
                        GCellGraph& graph,
                        std::string& error)
{
  std::string content;
  if (!readWholeFile(path, content, error)) {
    return false;
  }
  Cursor cursor(content);
  int64_t layers, x_size, y_size;
  if (!cursor.readInt64(layers) || !cursor.readInt64(x_size)
      || !cursor.readInt64(y_size)) {
    error = "bad resource header";
    return false;
  }
  graph = GCellGraph();
  graph.num_layers = static_cast<int>(layers);
  graph.x_size = static_cast<int>(x_size);
  graph.y_size = static_cast<int>(y_size);
  graph.min_routing_layer = 1;  // layer 0 holds pins only
  graph.layers.resize(layers);
  if (!cursor.readDouble(graph.unit_wire_cost)
      || !cursor.readDouble(graph.unit_via_cost)) {
    error = "bad unit costs";
    return false;
  }
  for (int z = 0; z < layers; ++z) {
    if (!cursor.readDouble(graph.layers[z].overflow_weight)) {
      error = "bad overflow weights";
      return false;
    }
  }
  graph.horizontal_edge_length.resize(x_size - 1);
  for (int64_t& length : graph.horizontal_edge_length) {
    if (!cursor.readInt64(length)) {
      error = "bad horizontal edge lengths";
      return false;
    }
  }
  graph.vertical_edge_length.resize(y_size - 1);
  for (int64_t& length : graph.vertical_edge_length) {
    if (!cursor.readInt64(length)) {
      error = "bad vertical edge lengths";
      return false;
    }
  }
  cursor.skipLineRemainder();
  graph.capacity.assign(layers, std::vector<float>(x_size * y_size, 0.0f));
  for (int z = 0; z < layers; ++z) {
    cursor.skipSpace();
    const std::string header = cursor.readLine();
    char name[256] = {0};
    int direction = 0;
    double min_length = 0;
    if (std::sscanf(
            header.c_str(), "%255s %d %lf", name, &direction, &min_length)
        != 3) {
      error = "bad layer header: " + header;
      return false;
    }
    RoutingLayer& layer = graph.layers[z];
    layer.name = name;
    layer.horizontal = direction == 0;
    layer.min_length = min_length;
    std::vector<float>& capacity = graph.capacity[z];
    for (int y = 0; y < y_size; ++y) {
      for (int x = 0; x < x_size; ++x) {
        double value;
        if (!cursor.readDouble(value)) {
          error = "bad capacity on layer " + layer.name;
          return false;
        }
        capacity[graph.layerEdgeIndex(z, x, y)] = static_cast<float>(value);
      }
    }
    cursor.skipLineRemainder();
  }
  return true;
}

bool readIspd24Nets(const std::string& path,
                    std::vector<RoutingNet>& nets,
                    std::string& error)
{
  std::string content;
  if (!readWholeFile(path, content, error)) {
    return false;
  }
  nets.clear();
  const char* p = content.data();
  const char* end = p + content.size();
  RoutingNet* net = nullptr;
  while (p < end) {
    const char* line_end
        = static_cast<const char*>(std::memchr(p, '\n', end - p));
    if (line_end == nullptr) {
      line_end = end;
    }
    const char* q = p;
    while (q < line_end && (*q == ' ' || *q == '\t' || *q == '\r')) {
      ++q;
    }
    if (q < line_end) {
      if (*q == '[') {
        if (net == nullptr) {
          error = "pin outside a net";
          return false;
        }
        RoutingPin& pin = net->pins.emplace_back();
        // Triples (layer, x, y) separated by punctuation.
        int values[3];
        int count = 0;
        while (q < line_end) {
          if (*q >= '0' && *q <= '9') {
            int value = 0;
            while (q < line_end && *q >= '0' && *q <= '9') {
              value = value * 10 + (*q - '0');
              ++q;
            }
            values[count++] = value;
            if (count == 3) {
              pin.access_points.push_back({values[0], values[1], values[2]});
              count = 0;
            }
          } else {
            ++q;
          }
        }
      } else if (*q == '(' || *q == ')') {
        if (*q == ')') {
          net = nullptr;
        }
      } else {
        const char* name_end = line_end;
        while (name_end > q
               && (name_end[-1] == ' ' || name_end[-1] == '\r'
                   || name_end[-1] == '\t')) {
          --name_end;
        }
        net = &nets.emplace_back();
        net->name.assign(q, name_end);
      }
    }
    p = line_end + 1;
  }
  return true;
}

bool writeIspd24Solution(const std::string& path,
                         const std::vector<RoutingNet>& nets,
                         const std::vector<NetRoute3D>& routes,
                         std::string& error)
{
  FILE* file = std::fopen(path.c_str(), "wb");
  if (file == nullptr) {
    error = "cannot open " + path + ": " + std::strerror(errno);
    return false;
  }
  std::vector<char> buffer(1 << 22);
  std::setvbuf(file, buffer.data(), _IOFBF, buffer.size());
  for (size_t n = 0; n < nets.size(); ++n) {
    std::fputs(nets[n].name.c_str(), file);
    std::fputs("\n(\n", file);
    for (const RouteSegment3D& s : routes[n].segments) {
      std::fprintf(file,
                   "%d %d %d %d %d %d\n",
                   s.x_low,
                   s.y_low,
                   s.layer_low,
                   s.x_high,
                   s.y_high,
                   s.layer_high);
    }
    std::fputs(")\n", file);
  }
  if (std::fclose(file) != 0) {
    error = "write failed on " + path;
    return false;
  }
  return true;
}

}  // namespace grt::helem
