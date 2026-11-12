// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

// Differential test of the BonnRoute routing-space rule model against drt's
// GC on asap7 rules (test/asap7/asap7_tech_1x_201209.lef, DBU 4000):
// LEF58_EOLKEEPOUT (with CORNERONLY), LEF58_CORNERSPACING, LEF58_SPACING
// ENDOFLINE ... ENDTOEND (with PARALLELEDGE on M2/M3) and LEF58_CUTCLASS +
// LEF58_SPACINGTABLE. Random pairs of BonnRoute's shape kinds (routed wire
// with its via pads, via pad, pin, obstruction, cut), of different nets
// and of one net, are placed far apart in one design; drt's markers per
// pair are compared with BonnRoute's verdict on the same geometry
// (br::violationMask with the shape classes BonnRouteDR uses, rules from
// br::drtWiringRules / drtViaRules). A pair where BonnRoute says legal and
// drt reports a violation is an optimistic mismatch (a failure); the
// reverse is conservative and only counted. Tables go to stdout.

#include <algorithm>
#include <array>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "br/BrDrtRules.h"
#include "br/BrGeometry.h"
#include "br/BrWireModel.h"
#include "db/obj/frMarker.h"
#include "db/obj/frNet.h"
#include "db/obj/frVia.h"
#include "db/tech/frConstraint.h"
#include "db/tech/frLayer.h"
#include "db/tech/frLookupTbl.h"
#include "db/tech/frViaDef.h"
#include "fixture.h"
#include "frBaseTypes.h"
#include "frDesign.h"
#include "gc/FlexGC.h"
#include "gtest/gtest.h"
#include "odb/db.h"
#include "odb/dbTypes.h"
#include "odb/geom.h"

namespace drt {
namespace {

using odb::dbTechLayerDir;
using odb::dbTechLayerType;

constexpr frCoord kCell = 8000;  // pairs are this far apart

enum class Kind : uint8_t
{
  kWire,  // routed preferred-direction wire with a via pad at each end
  kPad,   // via pad
  kPin,   // instance pin (fixed)
  kObs,   // obstruction (fixed)
  kCut,   // via (cut with its pads; only cut-layer markers count)
};
const char* kindName(const Kind k)
{
  switch (k) {
    case Kind::kWire:
      return "wire+pads";
    case Kind::kPad:
      return "pad";
    case Kind::kPin:
      return "pin";
    case Kind::kObs:
      return "obs";
    case Kind::kCut:
      return "cut";
  }
  return "?";
}

struct Shape
{
  frLayerNum layer;
  odb::Rect rect;
  br::ShapeRole role;
};

struct Item
{
  Kind kind = Kind::kPad;
  std::vector<Shape> shapes;  // relative to the item origin
  const frViaDef* via = nullptr;
  odb::Point via_at;  // relative
};

// drt marker families.
enum DrtFamily : uint32_t
{
  kDrtShort = 1u << 0,
  kDrtSpacing = 1u << 1,
  kDrtEol = 1u << 2,
  kDrtKeepOut = 1u << 3,
  kDrtCorner = 1u << 4,
  kDrtCutTable = 1u << 5,
  kDrtCut = 1u << 6,
  kDrtOther = 1u << 7,
};
constexpr int kDrtFamilies = 8;
const char* kDrtNames[kDrtFamilies] = {"Short",
                                       "MetalSpacing",
                                       "Lef58EOL",
                                       "EolKeepOut",
                                       "CornerSpacing",
                                       "CutSpcTbl",
                                       "CutSpacing",
                                       "other"};
constexpr uint32_t kDrtPairRules = ~kDrtOther;

uint32_t drtFamily(const frConstraint* con)
{
  if (con == nullptr) {
    return kDrtOther;
  }
  switch (con->typeId()) {
    case frConstraintTypeEnum::frcShortConstraint:
      return kDrtShort;
    case frConstraintTypeEnum::frcSpacingConstraint:
    case frConstraintTypeEnum::frcSpacingTablePrlConstraint:
    case frConstraintTypeEnum::frcSpacingTableTwConstraint:
    case frConstraintTypeEnum::frcSpacingRangeConstraint:
      return kDrtSpacing;
    case frConstraintTypeEnum::frcLef58SpacingEndOfLineConstraint:
      return kDrtEol;
    case frConstraintTypeEnum::frcLef58EolKeepOutConstraint:
      return kDrtKeepOut;
    case frConstraintTypeEnum::frcLef58CornerSpacingConstraint:
      return kDrtCorner;
    case frConstraintTypeEnum::frcLef58CutSpacingTableConstraint:
      return kDrtCutTable;
    case frConstraintTypeEnum::frcCutSpacingConstraint:
      return kDrtCut;
    default:
      return kDrtOther;
  }
}

// BonnRoute rule families: violationMask bits (spacing, line-end, the
// side rules and corner spacing when present) and the cut rule.
constexpr int kBrFamilies = 8;
const char* kBrNames[kBrFamilies] = {"spacing",
                                     "line-end",
                                     "keepout-region",
                                     "keepout-corner",
                                     "eol-facing",
                                     "end-to-end",
                                     "corner",
                                     "?"};

struct Pair
{
  int cell = 0;
  frLayerNum layer = 0;  // layer whose markers count
  bool same_net = false;
  bool control = false;  // item a alone (self-violation control)
  Item a;
  Item b;
  odb::Point oa;  // absolute origins
  odb::Point ob;
  uint32_t drt = 0;
  uint32_t drt_other = 0;
  uint32_t br = 0;
};

struct LayerSpec
{
  const char* name;
  bool horizontal;
  frCoord width;
  frCoord pitch;
};

// asap7 (DBU 4000): M1..M6 and V1..V5.
constexpr LayerSpec kMetal[] = {{"m1", false, 72, 144},
                                {"m2", true, 72, 144},
                                {"m3", false, 72, 144},
                                {"m4", true, 96, 192},
                                {"m5", false, 96, 192},
                                {"m6", true, 128, 256}};

struct Lef58DiffFixture : public Fixture
{
  Lef58DiffFixture() { buildTech(); }

  using Fixture::db_tech;
  using Fixture::design;
  using Fixture::initRegionQuery;
  using Fixture::logger;
  using Fixture::makeInst;
  using Fixture::makeMacro;
  using Fixture::makeNet;
  using Fixture::router_cfg;

  frLayer* layer(const frLayerNum l) { return design->getTech()->getLayer(l); }

  void setWidth(frLayer* l, const frCoord w)
  {
    l->setWidth(w);
    l->setMinWidth(w);
    auto con = std::make_unique<frMinWidthConstraint>(w);
    l->setMinWidthConstraint(con.get());
    design->getTech()->addUConstraint(std::move(con));
  }

  void addMetal(const LayerSpec& s)
  {
    addLayer(
        design->getTech(),
        s.name,
        dbTechLayerType::ROUTING,
        s.horizontal ? dbTechLayerDir::HORIZONTAL : dbTechLayerDir::VERTICAL);
  }

  void spacingTable(const frLayerNum l,
                    const frCoord base,
                    const frCoord wide_width,
                    const frCoord wide)
  {
    fr2DLookupTbl<frCoord, frCoord, frCoord> tbl(
        "WIDTH", {0, wide_width}, "PARALLELRUNLENGTH", {0}, {{base}, {wide}});
    auto con = std::make_unique<frSpacingTablePrlConstraint>(tbl);
    layer(l)->setMinSpacing(con.get());
    design->getTech()->addUConstraint(std::move(con));
  }

  void eol(const frLayerNum l,
           const frCoord space,
           const frCoord width,
           const frCoord within,
           const frCoord ete,
           const frCoord par_space,
           const frCoord par_within)
  {
    frLef58SpacingEndOfLineConstraint* con
        = makeLef58SpacingEolConstraint(l, space, width, within);
    auto e
        = std::make_shared<frLef58SpacingEndOfLineWithinEndToEndConstraint>();
    e->setEndToEndSpace(ete);
    con->getWithinConstraint()->setEndToEndConstraint(e);
    if (par_space > 0) {
      makeLef58SpacingEolParEdgeConstraint(con, par_space, par_within);
    }
  }

  void corner(const frLayerNum l,
              const frCoord width,
              const frCoord spacing,
              const frCoord within)
  {
    fr1DLookupTbl<frCoord, std::pair<frCoord, frCoord>> tbl(
        "WIDTH", {width}, {{spacing, spacing}});
    auto con = std::make_unique<frLef58CornerSpacingConstraint>(tbl);
    con->setCornerType(frCornerTypeEnum::CONVEX);
    con->setSameXY(true);
    if (within >= 0) {
      con->setWithin(within);
    }
    layer(l)->addLef58CornerSpacingConstraint(con.get());
    design->getTech()->addUConstraint(std::move(con));
  }

  void cutTable(const frLayerNum l,
                const std::vector<std::pair<std::string, odb::Point>>& cls,
                const frCoord dflt,
                const std::pair<int, int> special,
                const frCoord special_spacing)
  {
    for (const auto& [name, wl] : cls) {
      makeCutClass(l, name, wl.x(), wl.y());
    }
    auto dbRule = odb::dbTechLayerCutSpacingTableDefRule::create(
        db_tech->findLayer(layer(l)->getName().c_str()));
    dbRule->setDefault(dflt);
    dbRule->setDefaultValid(true);
    // As the LEF parser stores the table: SIDE and END columns and rows
    // per class.
    const int n = static_cast<int>(cls.size());
    std::map<std::string, uint32_t> rows;
    std::map<std::string, uint32_t> cols;
    for (int k = 0; k < n; k++) {
      rows[cls[k].first + "/SIDE"] = k;
      rows[cls[k].first + "/END"] = n + k;
      cols[cls[k].first + "/SIDE"] = k;
      cols[cls[k].first + "/END"] = n + k;
    }
    std::vector<std::vector<std::pair<int, int>>> table(
        2 * n, std::vector<std::pair<int, int>>(2 * n, {dflt, dflt}));
    if (special.first >= 0) {
      for (const int r : {special.first, special.first + n}) {
        for (const int c : {special.second, special.second + n}) {
          table[r][c] = {special_spacing, special_spacing};
          table[c][r] = {special_spacing, special_spacing};
        }
      }
    }
    dbRule->setSpacingTable(table, rows, cols);
    makeLef58CutSpcTbl(l, dbRule);
  }

  frViaDef* via(const char* name,
                const frLayerNum cut,
                const odb::Rect& bottom,
                const odb::Rect& cut_rect,
                const odb::Rect& top)
  {
    auto v = std::make_unique<frViaDef>(name);
    auto fig = [&](const frLayerNum l, const odb::Rect& r) {
      auto f = std::make_unique<frRect>();
      f->setBBox(r);
      f->setLayerNum(l);
      return f;
    };
    v->addLayer1Fig(fig(cut - 1, bottom));
    v->addCutFig(fig(cut, cut_rect));
    v->addLayer2Fig(fig(cut + 1, top));
    return design->getTech()->addVia(std::move(v));
  }

  // drt's summary of the EOL rules as FlexRP::prep_eolForbiddenLen fills
  // it (the line-end rule of the routing space before the LEF58 model).
  void drEolSummary(frLayer* l)
  {
    frCoord eol = INT_MAX;
    for (const auto con : l->getLef58SpacingEndOfLineConstraints()) {
      eol = std::min(eol, con->getEolWidth());
    }
    for (const auto con : l->getLef58EolKeepOutConstraints()) {
      eol = std::min(eol, con->getEolWidth());
    }
    const frCoord width = eol == INT_MAX
                              ? static_cast<frCoord>(l->getWidth())
                              : std::max<frCoord>(eol - 1, l->getWidth());
    frCoord space = 0;
    frCoord within = 0;
    for (const auto con : l->getLef58SpacingEndOfLineConstraints()) {
      if (width < con->getEolWidth()) {
        space = std::max(space, con->getEolSpace());
        if (con->hasWithinConstraint()) {
          const auto w = con->getWithinConstraint();
          within = std::max(within, w->getEolWithin());
          if (w->hasEndToEndConstraint()) {
            space = std::max(space,
                             w->getEndToEndConstraint()->getEndToEndSpace());
          }
        }
      }
    }
    for (const auto con : l->getLef58EolKeepOutConstraints()) {
      if (width < con->getEolWidth()) {
        space = std::max(space, con->getForwardExt());
        within = std::max(within, con->getSideExt());
      }
    }
    l->setDrEolSpacingConstraint(width, space, within);
  }

  void buildTech()
  {
    frTechObject* tech = design->getTech();
    tech->setManufacturingGrid(1);
    // The fixture made masterslice (0), v0 (1), m1 (2).
    layer(2)->getDbLayer()->setDirection(dbTechLayerDir::VERTICAL);
    for (int k = 1; k < 6; k++) {
      const std::string cut = "v" + std::to_string(k);
      addLayer(tech, cut.c_str(), dbTechLayerType::CUT);
      addMetal(kMetal[k]);
    }
    for (int k = 0; k < 6; k++) {
      frLayer* l = layer(2 + 2 * k);
      setWidth(l, kMetal[k].width);
      l->getDbLayer()->setPitch(kMetal[k].pitch);
    }
    // M1: SPACING 0.018; EOLKEEPOUT 0.01825 EXTENSION 0 0 0.031;
    // CORNERSPACING CONVEXCORNER CORNERONLY 0.01 WIDTH 0.018 SPACING 0.018.
    makeSimpleSpacingConstraint(2, 72);
    makeLef58EolKeepOutConstraint(2, false, false, -10, 10, 124, 0, 0, 73);
    corner(2, 72, 72, 40);
    // M2 / M3: SPACING 0.018 ENDOFLINE 0.025 WITHIN 0.02 (M3 0.0125)
    // ENDTOEND 0.031 PARALLELEDGE 0.025 WITHIN 0.02; EOLKEEPOUT 0.025
    // EXTENSION 0 0.0125 0.031 CORNERONLY; CORNERSPACING WIDTH 0 SPACING
    // 0.02.
    for (const frLayerNum l : {4, 6}) {
      makeSimpleSpacingConstraint(l, 72);
      eol(l, 72, 100, l == 4 ? 80 : 50, 124, 100, 80);
      makeLef58EolKeepOutConstraint(l, true, false, -10, 10, 124, 50, 0, 100);
      corner(l, 0, 80, -1);
    }
    // M4 / M5: SPACINGTABLE 0.024 / 0.072 above width 0.025; SPACING 0.024
    // ENDOFLINE 0.025 WITHIN 0.04 ENDTOEND 0.04; EOLKEEPOUT 0.025
    // EXTENSION 0.048 0.02425 0.048 CORNERONLY; CORNERSPACING CORNERONLY
    // 0.048 WIDTH 0 SPACING 0.04.
    for (const frLayerNum l : {8, 10}) {
      spacingTable(l, 96, 100, 288);
      eol(l, 96, 100, 160, 160, 0, 0);
      makeLef58EolKeepOutConstraint(l, true, false, -10, 10, 192, 97, 192, 100);
      corner(l, 0, 160, 192);
    }
    // M6: 0.032 / 0.072 above 0.033; EOL 0.032 0.0375 WITHIN 0.04
    // ENDTOEND 0.04; EOLKEEPOUT 0.05 EXTENSION 0.048 0.03225 0.048
    // CORNERONLY; CORNERSPACING as M4.
    spacingTable(12, 128, 132, 288);
    eol(12, 128, 150, 160, 160, 0, 0);
    makeLef58EolKeepOutConstraint(12, true, false, -10, 10, 192, 129, 192, 200);
    corner(12, 0, 160, 192);
    // V3: CUTCLASS V3 0.018 x 0.024, V3_0p480 0.018 x 0.12, V3_0p864
    // 0.018 x 0.216; SPACINGTABLE DEFAULT 0.034 (all entries "-").
    cutTable(
        7,
        {{"V3", {72, 96}}, {"V3_0p480", {72, 480}}, {"V3_0p864", {72, 864}}},
        136,
        {-1, -1},
        0);
    cutTable(
        9,
        {{"Vx", {96, 96}}, {"Vx_0p480", {96, 480}}, {"Vx_0p864", {96, 864}}},
        136,
        {-1, -1},
        0);
    // V5 as asap7 plus one class-pair entry (class lookup coverage).
    cutTable(
        11,
        {{"Vx", {96, 128}}, {"Vx_0p480", {96, 640}}, {"Vx_0p864", {96, 1152}}},
        136,
        {0, 1},
        200);
    // Default vias (asap7 VIA12 .. VIA56).
    vias_[3] = via(
        "VIA12", 3, {-36, -44, 36, 44}, {-36, -36, 36, 36}, {-56, -36, 56, 36});
    vias_[5] = via(
        "VIA23", 5, {-56, -36, 56, 36}, {-36, -36, 36, 36}, {-36, -56, 36, 56});
    vias_[7] = via(
        "VIA34", 7, {-36, -68, 36, 68}, {-36, -48, 36, 48}, {-80, -48, 80, 48});
    vias_[9] = via(
        "VIA45", 9, {-92, -48, 92, 48}, {-48, -48, 48, 48}, {-48, -92, 48, 92});
    vias_[11] = via("VIA56",
                    11,
                    {-48, -108, 48, 108},
                    {-48, -64, 48, 64},
                    {-92, -64, 92, 64});
    for (int k = 0; k < 6; k++) {
      drEolSummary(layer(2 + 2 * k));
    }
  }

  // BonnRoute rules of fr layers 2 (m1) .. 12 (m6) at index fr - 2.
  void buildRoutingSpace()
  {
    std::vector<br::LayerRules> rules;
    for (frLayerNum l = 2; l <= 12; l++) {
      frLayer* fl = layer(l);
      if (fl->getType() == dbTechLayerType::ROUTING) {
        rules.push_back(br::drtWiringRules(fl));
      } else {
        const frVia v(vias_[l]);
        frCoord s = fl->hasCutSpacing() ? fl->getCutSpacingValue() : 0;
        if (s <= 0) {
          s = v.getCutBBox().minDXDY();
        }
        rules.push_back(
            br::drtViaRules(fl, nullptr, s, layer(l - 1)->isHorizontal()));
      }
    }
    classes_.setLayers(rules);
  }

  br::ClassId classOf(const Shape& s)
  {
    const br::Rect r{
        s.rect.xMin(), s.rect.yMin(), s.rect.xMax(), s.rect.yMax()};
    return classes_.classify(s.layer - 2, r, s.role);
  }

  // Pads of the vias landing on metal layer l (relative to the via point).
  std::vector<odb::Rect> padsOn(const frLayerNum l)
  {
    std::vector<odb::Rect> out;
    if (vias_[l - 1] != nullptr) {
      out.push_back(frVia(vias_[l - 1]).getLayer2BBox());
    }
    if (l + 1 <= 11 && vias_[l + 1] != nullptr) {
      out.push_back(frVia(vias_[l + 1]).getLayer1BBox());
    }
    return out;
  }

  static int rnd(std::mt19937& rng, const int lo, const int hi)
  {
    return lo + static_cast<int>(rng() % static_cast<unsigned>(hi - lo + 1));
  }

  Item makeItem(std::mt19937& rng, const frLayerNum l, const Kind kind)
  {
    Item it;
    it.kind = kind;
    const LayerSpec& s = kMetal[(l - 2) / 2];
    const bool h = s.horizontal;
    auto oriented = [&](const frCoord a_lo,
                        const frCoord a_hi,
                        const frCoord c_lo,
                        const frCoord c_hi) {
      return h ? odb::Rect(a_lo, c_lo, a_hi, c_hi)
               : odb::Rect(c_lo, a_lo, c_hi, a_hi);
    };
    const std::vector<odb::Rect> pads = padsOn(l);
    auto pad_at = [&](const frCoord along) {
      odb::Rect p = pads[rng() % pads.size()];
      p.moveDelta(h ? along : 0, h ? 0 : along);
      return p;
    };
    switch (kind) {
      case Kind::kWire: {
        const frCoord w = s.width;
        const frCoord len = rnd(rng, 1, 8) * s.pitch + rnd(rng, 0, 1) * 48;
        it.shapes.push_back(
            {l,
             oriented(-w / 2, len + w / 2, -w / 2, w / 2),
             l == 2 ? br::ShapeRole::kFixed : br::ShapeRole::kPreferredWire});
        it.shapes.push_back({l, pad_at(0), br::ShapeRole::kViaPad});
        it.shapes.push_back({l, pad_at(len), br::ShapeRole::kViaPad});
        break;
      }
      case Kind::kPad:
        it.shapes.push_back({l, pad_at(0), br::ShapeRole::kViaPad});
        break;
      case Kind::kPin: {
        const frCoord w = rng() % 3 == 0 ? rnd(rng, 36, 160) : s.width;
        const frCoord len = rnd(rng, 80, 600);
        it.shapes.push_back(
            {l, oriented(0, len, -w / 2, w - w / 2), br::ShapeRole::kFixed});
        break;
      }
      case Kind::kObs:
      case Kind::kCut: {
        const frCoord dx = rnd(rng, 36, 600);
        const frCoord dy = rnd(rng, 36, 600);
        it.shapes.push_back(
            {l, odb::Rect(0, 0, dx, dy), br::ShapeRole::kFixed});
        break;
      }
    }
    return it;
  }

  Item makeCutItem(std::mt19937& rng, const frLayerNum cut)
  {
    Item it;
    it.kind = Kind::kCut;
    const frViaDef* def = vias_[cut];
    const int pick = rng() % 8;
    if (pick < 5) {
      it.via = def;
    } else {
      // Another cut class (or none): a cut of that size with the default
      // via's pads.
      static const frCoord kLen[] = {480, 864, 72};
      const frCoord w = frVia(def).getCutBBox().minDXDY();
      const frCoord len = kLen[pick - 5];
      const std::string name
          = "cut" + std::to_string(cut) + "_" + std::to_string(len);
      auto key = std::make_pair(cut, len);
      if (!extra_vias_.contains(key)) {
        const frVia v(def);
        const odb::Rect c(-w / 2, -len / 2, w - w / 2, len - len / 2);
        odb::Rect b = v.getLayer1BBox();
        odb::Rect t = v.getLayer2BBox();
        b.merge(c);
        t.merge(c);
        extra_vias_[key] = via(name.c_str(), cut, b, c, t);
      }
      it.via = extra_vias_[key];
    }
    it.shapes.push_back({cut, frVia(it.via).getCutBBox(), br::ShapeRole::kCut});
    return it;
  }

  static bool touches(const Item& a,
                      const odb::Point oa,
                      const Item& b,
                      const odb::Point ob)
  {
    for (const Shape& s : a.shapes) {
      odb::Rect ra = s.rect;
      ra.moveDelta(oa.x(), oa.y());
      for (const Shape& t : b.shapes) {
        if (t.layer != s.layer) {
          continue;
        }
        odb::Rect rb = t.rect;
        rb.moveDelta(ob.x(), ob.y());
        if (ra.intersects(rb)) {
          return true;
        }
      }
    }
    return false;
  }

  void addItem(const Item& it,
               const odb::Point o,
               frNet* net,
               frMaster* master,
               std::vector<std::pair<int, frNet*>>& pin_nets)
  {
    if (it.via != nullptr) {
      makeVia(const_cast<frViaDef*>(it.via), net, o);
      return;
    }
    for (const Shape& s : it.shapes) {
      odb::Rect r = s.rect;
      r.moveDelta(o.x(), o.y());
      if (it.kind == Kind::kPin) {
        makeMacroPin(master,
                     "p" + std::to_string(pin_nets.size()),
                     r.xMin(),
                     r.yMin(),
                     r.xMax(),
                     r.yMax(),
                     s.layer);
        pin_nets.emplace_back(static_cast<int>(pin_nets.size()), net);
        continue;
      }
      if (it.kind == Kind::kObs) {
        makeMacroObs(master, r.xMin(), r.yMin(), r.xMax(), r.yMax(), s.layer);
        continue;
      }
      // Metal of a routed net: a path segment with truncated ends.
      const bool hseg = r.dx() >= r.dy();
      const frCoord w = hseg ? r.dy() : r.dx();
      const odb::Point c = r.center();
      if (hseg) {
        makePathseg(net, s.layer, {r.xMin(), c.y()}, {r.xMax(), c.y()}, w);
      } else {
        makePathseg(net, s.layer, {c.x(), r.yMin()}, {c.x(), r.yMax()}, w);
      }
    }
  }

  uint32_t brVerdict(const Pair& p)
  {
    uint32_t mask = 0;
    for (const Shape& sa : p.a.shapes) {
      if (sa.layer != p.layer) {
        continue;
      }
      const br::Rect ra = brRect(sa.rect, p.oa);
      const br::ClassId ca = classOf(sa);
      for (const Shape& sb : p.b.shapes) {
        if (sb.layer != p.layer) {
          continue;
        }
        mask |= br::violationMask(
            classes_,
            br::PlacedShape{ra, ca, 3},
            br::PlacedShape{brRect(sb.rect, p.ob), classOf(sb), 3});
      }
    }
    return mask;
  }

  static br::Rect brRect(const odb::Rect& r, const odb::Point o)
  {
    return {
        r.xMin() + o.x(), r.yMin() + o.y(), r.xMax() + o.x(), r.yMax() + o.y()};
  }

  std::array<frViaDef*, 13> vias_{};
  std::map<std::pair<frLayerNum, frCoord>, frViaDef*> extra_vias_;
  br::ShapeClassTable classes_;
};

std::string layerName(const frLayerNum l)
{
  return (l % 2 == 0 ? "M" : "V") + std::to_string(l / 2);
}

TEST_F(Lef58DiffFixture, asap7_pairs)
{
  Lef58DiffFixture& f = *this;
  f.buildRoutingSpace();
  const char* env = std::getenv("BR_DIFF_SEED");
  std::mt19937 rng(env != nullptr ? std::atoi(env) : 20261112);
  const char* per_env = std::getenv("BR_DIFF_PAIRS");
  const int per = per_env != nullptr ? std::atoi(per_env) : 120;

  std::vector<Pair> pairs;
  const Kind routed[] = {Kind::kWire, Kind::kPad};
  const Kind others[] = {Kind::kWire, Kind::kPad, Kind::kPin, Kind::kObs};
  for (frLayerNum l = 2; l <= 12; l += 2) {
    const LayerSpec& s = kMetal[(l - 2) / 2];
    for (const Kind ka : routed) {
      for (const Kind kb : others) {
        for (const bool same : {false, true}) {
          if (same && kb == Kind::kObs) {
            continue;
          }
          for (int k = 0; k < per; k++) {
            Pair p;
            p.layer = l;
            p.same_net = same;
            p.a = f.makeItem(rng, l, ka);
            p.b = f.makeItem(rng, l, kb);
            bool placed = false;
            for (int attempt = 0; attempt < 50 && !placed; attempt++) {
              frCoord dx;
              frCoord dy;
              if (rng() % 2 == 0) {
                // Track-aligned: 0..2 tracks across, random along.
                const frCoord across
                    = Lef58DiffFixture::rnd(rng, -2, 2) * s.pitch;
                const frCoord along = Lef58DiffFixture::rnd(rng, -900, 900);
                dx = s.horizontal ? along : across;
                dy = s.horizontal ? across : along;
              } else {
                dx = Lef58DiffFixture::rnd(rng, -700, 700);
                dy = Lef58DiffFixture::rnd(rng, -700, 700);
              }
              p.ob = odb::Point(dx, dy);
              placed
                  = !same || !Lef58DiffFixture::touches(p.a, {0, 0}, p.b, p.ob);
            }
            if (placed) {
              pairs.push_back(std::move(p));
            }
          }
        }
      }
    }
    // Controls: item a alone.
    for (const Kind ka : routed) {
      for (int k = 0; k < per / 4; k++) {
        Pair p;
        p.layer = l;
        p.control = true;
        p.a = f.makeItem(rng, l, ka);
        pairs.push_back(std::move(p));
      }
    }
  }
  for (const frLayerNum cut : {7, 9, 11}) {
    const LayerSpec& below = kMetal[(cut - 3) / 2];
    const LayerSpec& above = kMetal[(cut - 1) / 2];
    const frCoord px = below.horizontal ? above.pitch : below.pitch;
    const frCoord py = below.horizontal ? below.pitch : above.pitch;
    for (const bool same : {false, true}) {
      for (int k = 0; k < 4 * per; k++) {
        Pair p;
        p.layer = cut;
        p.same_net = same;
        p.a = f.makeCutItem(rng, cut);
        p.b = f.makeCutItem(rng, cut);
        frCoord dx = Lef58DiffFixture::rnd(rng, -2, 2) * px;
        frCoord dy = Lef58DiffFixture::rnd(rng, -2, 2) * py;
        if (rng() % 3 == 0) {
          dx += Lef58DiffFixture::rnd(rng, -120, 120);
          dy += Lef58DiffFixture::rnd(rng, -120, 120);
        }
        if (dx == 0 && dy == 0) {
          dx = px;
        }
        p.ob = odb::Point(dx, dy);
        pairs.push_back(std::move(p));
      }
    }
  }

  // One design: pair k in cell k.
  const int cols = 160;
  auto [master, db_master] = f.makeMacro("PINS", 0, 0, kCell, kCell);
  std::vector<std::pair<int, frNet*>> pin_nets;
  for (size_t k = 0; k < pairs.size(); k++) {
    Pair& p = pairs[k];
    p.cell = static_cast<int>(k);
    const int cell = static_cast<int>(k);
    const odb::Point o((cell % cols) * kCell + kCell / 2,
                       (cell / cols) * kCell + kCell / 2);
    p.oa = o;
    p.ob = odb::Point(o.x() + p.ob.x(), o.y() + p.ob.y());
    frNet* na = f.makeNet(("a" + std::to_string(k)).c_str());
    f.addItem(p.a, p.oa, na, master, pin_nets);
    if (!p.control) {
      frNet* nb
          = p.same_net ? na : f.makeNet(("b" + std::to_string(k)).c_str());
      f.addItem(p.b, p.ob, nb, master, pin_nets);
    }
  }
  frInst* inst = f.makeInst("pins", master, db_master);
  for (const auto& [term, net] : pin_nets) {
    frInstTerm* it = inst->getInstTerms()[term].get();
    it->addToNet(net);
    net->addInstTerm(it);
  }
  f.initRegionQuery();
  const int rows = static_cast<int>((pairs.size() + cols - 1) / cols);
  const odb::Rect work(0, 0, cols * kCell, rows * kCell);
  FlexGCWorker worker(f.design->getTech(), f.logger.get(), f.router_cfg.get());
  worker.setExtBox(work);
  worker.setDrcBox(work);
  worker.init(f.design.get());
  worker.main();
  int64_t stray = 0;
  for (const auto& m : worker.getMarkers()) {
    const odb::Rect b = m->getBBox();
    const int cx = (b.xMin() + b.xMax()) / 2 / kCell;
    const int cy = (b.yMin() + b.yMax()) / 2 / kCell;
    const size_t k = static_cast<size_t>(cy) * cols + cx;
    if (k >= pairs.size()) {
      stray++;
      continue;
    }
    Pair& p = pairs[k];
    if (m->getLayerNum() != p.layer) {
      continue;
    }
    const uint32_t fam = drtFamily(m->getConstraint());
    if (fam == kDrtOther) {
      p.drt_other |= 1;
    } else {
      p.drt |= fam;
    }
  }
  for (Pair& p : pairs) {
    if (!p.control) {
      p.br = f.brVerdict(p);
    }
  }

  // Tallies.
  struct Row
  {
    int64_t pairs = 0, drt = 0, br = 0, both = 0, optimistic = 0,
            conservative = 0, other = 0;
  };
  std::map<std::pair<frLayerNum, int>, Row> by_layer;
  std::map<std::tuple<frLayerNum, int, int>, int64_t> opt_by_rule;
  std::map<std::tuple<frLayerNum, int, int>, int64_t> cons_by_rule;
  std::map<std::tuple<frLayerNum, std::string, int>, int64_t> opt_by_kind;
  int64_t control_markers = 0;
  int64_t optimistic = 0;
  int printed = 0;
  for (const Pair& p : pairs) {
    if (p.control) {
      control_markers += (p.drt != 0) + p.drt_other;
      continue;
    }
    Row& r = by_layer[{p.layer, p.same_net}];
    r.pairs++;
    r.other += p.drt_other;
    const bool d = (p.drt & kDrtPairRules) != 0;
    const bool b = p.br != 0;
    r.drt += d;
    r.br += b;
    r.both += d && b;
    if (d && !b) {
      r.optimistic++;
      optimistic++;
      for (int k = 0; k < kDrtFamilies; k++) {
        if (p.drt & (1u << k)) {
          opt_by_rule[{p.layer, p.same_net, k}]++;
        }
      }
      opt_by_kind[{p.layer,
                   std::string(kindName(p.a.kind)) + " / " + kindName(p.b.kind),
                   p.same_net}]++;
      if (printed++ < 40) {
        std::printf("OPTIMISTIC %s %s drt=0x%x a=%s b=%s:",
                    layerName(p.layer).c_str(),
                    p.same_net ? "same-net" : "diff-net",
                    p.drt,
                    kindName(p.a.kind),
                    kindName(p.b.kind));
        for (const Item* it : {&p.a, &p.b}) {
          const odb::Point o = it == &p.a ? p.oa : p.ob;
          for (const Shape& s : it->shapes) {
            if (s.layer == p.layer) {
              std::printf(" [%d %d %d %d]",
                          s.rect.xMin() + o.x() - p.oa.x(),
                          s.rect.yMin() + o.y() - p.oa.y(),
                          s.rect.xMax() + o.x() - p.oa.x(),
                          s.rect.yMax() + o.y() - p.oa.y());
            }
          }
          std::printf(" |");
        }
        std::printf("\n");
      }
    }
    if (b && !d) {
      r.conservative++;
      for (int k = 0; k < kBrFamilies; k++) {
        if (p.br & (1u << k)) {
          cons_by_rule[{p.layer, p.same_net, k}]++;
        }
      }
    }
  }
  std::printf(
      "\n## Pairs by layer (drt = drt reports a pair rule; br = BonnRoute "
      "illegal)\n\n");
  std::printf(
      "| layer | nets | pairs | drt | br | both | optimistic | conservative | "
      "drt other markers |\n");
  std::printf("|---|---|---:|---:|---:|---:|---:|---:|---:|\n");
  Row total;
  for (const auto& [key, r] : by_layer) {
    std::printf(
        "| %s | %s | %lld | %lld | %lld | %lld | %lld | %lld | %lld |\n",
        layerName(key.first).c_str(),
        key.second ? "same" : "diff",
        (long long) r.pairs,
        (long long) r.drt,
        (long long) r.br,
        (long long) r.both,
        (long long) r.optimistic,
        (long long) r.conservative,
        (long long) r.other);
    total.pairs += r.pairs;
    total.drt += r.drt;
    total.br += r.br;
    total.both += r.both;
    total.optimistic += r.optimistic;
    total.conservative += r.conservative;
    total.other += r.other;
  }
  std::printf(
      "| all | all | %lld | %lld | %lld | %lld | %lld | %lld | %lld |\n",
      (long long) total.pairs,
      (long long) total.drt,
      (long long) total.br,
      (long long) total.both,
      (long long) total.optimistic,
      (long long) total.conservative,
      (long long) total.other);
  std::printf(
      "\n## Optimistic mismatches by drt rule (BonnRoute legal, drt flags the "
      "rule)\n\n");
  std::printf("| layer | nets | drt rule | pairs |\n|---|---|---|---:|\n");
  std::map<int, int64_t> opt_rule_total;
  for (const auto& [key, n] : opt_by_rule) {
    std::printf("| %s | %s | %s | %lld |\n",
                layerName(std::get<0>(key)).c_str(),
                std::get<1>(key) ? "same" : "diff",
                kDrtNames[std::get<2>(key)],
                (long long) n);
    opt_rule_total[std::get<2>(key)] += n;
  }
  for (const auto& [k, n] : opt_rule_total) {
    std::printf("| all | all | %s | %lld |\n", kDrtNames[k], (long long) n);
  }
  std::printf("\n## Optimistic mismatches by shape kinds\n\n");
  std::printf("| layer | kinds | nets | pairs |\n|---|---|---|---:|\n");
  for (const auto& [key, n] : opt_by_kind) {
    std::printf("| %s | %s | %s | %lld |\n",
                layerName(std::get<0>(key)).c_str(),
                std::get<1>(key).c_str(),
                std::get<2>(key) ? "same" : "diff",
                (long long) n);
  }
  std::printf(
      "\n## Conservative mismatches by BonnRoute rule (BonnRoute illegal, drt "
      "legal)\n\n");
  std::printf(
      "| layer | nets | BonnRoute rule | pairs |\n|---|---|---|---:|\n");
  std::map<int, int64_t> cons_rule_total;
  for (const auto& [key, n] : cons_by_rule) {
    std::printf("| %s | %s | %s | %lld |\n",
                layerName(std::get<0>(key)).c_str(),
                std::get<1>(key) ? "same" : "diff",
                kBrNames[std::get<2>(key)],
                (long long) n);
    cons_rule_total[std::get<2>(key)] += n;
  }
  for (const auto& [k, n] : cons_rule_total) {
    std::printf("| all | all | %s | %lld |\n", kBrNames[k], (long long) n);
  }
  std::printf("\ncontrol markers (item alone): %lld; stray markers: %lld\n",
              (long long) control_markers,
              (long long) stray);
  std::printf("TOTAL pairs=%lld optimistic=%lld conservative=%lld\n",
              (long long) total.pairs,
              (long long) total.optimistic,
              (long long) total.conservative);
  EXPECT_EQ(control_markers, 0);
  EXPECT_EQ(stray, 0);
  EXPECT_EQ(optimistic, 0);
}

}  // namespace
}  // namespace drt
