// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "br/BrIoPinRepair.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <set>
#include <utility>
#include <vector>

#include "br/BrIoPinArea.h"
#include "db/obj/frBTerm.h"
#include "db/obj/frBlock.h"
#include "db/obj/frNet.h"
#include "db/obj/frShape.h"
#include "db/obj/frVia.h"
#include "db/tech/frConstraint.h"
#include "db/tech/frLayer.h"
#include "db/tech/frTechObject.h"
#include "db/tech/frViaDef.h"
#include "drt-global.h"
#include "frBaseTypes.h"
#include "frDesign.h"
#include "frRegionQuery.h"
#include "gc/FlexGC.h"
#include "odb/dbTypes.h"
#include "odb/geom.h"
#include "utl/Logger.h"

namespace drt::br {
namespace {

Box toBox(const odb::Rect& box)
{
  return {box.xMin(), box.yMin(), box.xMax(), box.yMax()};
}

// Compare the real rule checker with and without a prospective patch.
// initPA0 reads the existing design; addPAObj adds only the candidate to
// the checker's temporary polygons. The routed design is unchanged.
bool patchClean(frDesign* design,
                 utl::Logger* logger,
                 RouterConfiguration* router_cfg,
                 frNet* net,
                 frPatchWire& patch)
{
  odb::Rect drc_box;
  odb::Rect ext_box;
  patch.getBBox().bloat(router_cfg->DRCSAFEDIST, drc_box);
  patch.getBBox().bloat(router_cfg->MTSAFEDIST, ext_box);
  std::set<MarkerId> before;
  for (const bool candidate : {false, true}) {
    FlexGCWorker checker(design->getTech(), logger, router_cfg);
    checker.setDrcBox(drc_box);
    checker.setExtBox(ext_box);
    checker.initPA0(design);
    if (candidate) {
      checker.addPAObj(&patch, net);
    }
    checker.initPA1();
    checker.setTargetNet(net);
    checker.main();
    for (const auto& marker : checker.getMarkers()) {
      const MarkerId key{marker->getBBox(),
                         marker->getLayerNum(),
                         marker->getConstraint(),
                         marker->getSrcs()};
      if (candidate && !before.contains(key)) {
        return false;
      }
      if (!candidate) {
        before.insert(key);
      }
    }
  }
  return true;
}

struct LayerMetal
{
  std::vector<Box> pieces;
  std::vector<odb::Point> anchors;
};

LayerMetal layerMetal(const frNet* net, const frLayerNum layer)
{
  LayerMetal result;
  // Surviving via origins first, then surviving wire ends.
  for (const auto& via : net->getVias()) {
    const frViaDef* def = via->getViaDef();
    if (def->getLayer1Num() == layer) {
      result.pieces.push_back(toBox(via->getLayer1BBox()));
      result.anchors.push_back(via->getOrigin());
    }
    if (def->getLayer2Num() == layer) {
      result.pieces.push_back(toBox(via->getLayer2BBox()));
      result.anchors.push_back(via->getOrigin());
    }
  }
  for (const auto& shape : net->getShapes()) {
    if (shape->typeId() != frcPathSeg || shape->getLayerNum() != layer) {
      continue;
    }
    const auto* seg = static_cast<const frPathSeg*>(shape.get());
    result.pieces.push_back(toBox(seg->getBBox()));
    const auto [begin, end] = seg->getPoints();
    result.anchors.push_back(begin);
    result.anchors.push_back(end);
  }
  for (const auto& patch : net->getPatchWires()) {
    if (patch->getLayerNum() == layer) {
      result.pieces.push_back(toBox(patch->getBBox()));
    }
  }
  return result;
}

}  // namespace

void repairIoPinArea(frDesign* design,
                     utl::Logger* logger,
                     RouterConfiguration* router_cfg)
{
  int64_t components = 0;
  int64_t below = 0;
  int64_t patched = 0;
  const GeoCoord grid
      = std::max<GeoCoord>(1, design->getTech()->getManufacturingGrid());
  const Box die = toBox(design->getTopBlock()->getDieBox());
  for (const auto& owned_net : design->getTopBlock()->getNets()) {
    frNet* net = owned_net.get();
    if (net->isSpecial() || net->isFake() || net->isFixed()
        || net->getType() == odb::dbSigType::POWER
        || net->getType() == odb::dbSigType::GROUND
        || net->getBTerms().empty()) {
      continue;
    }
    std::vector<frRect> pins;
    for (frBTerm* pin : net->getBTerms()) {
      pin->getShapes(pins);
    }
    for (const frRect& pin : pins) {
      const frLayerNum layer_num = pin.getLayerNum();
      frLayer* layer = design->getTech()->getLayer(layer_num);
      const frAreaConstraint* area = layer->getAreaConstraint();
      if (area == nullptr || area->getMinArea() <= 0) {
        continue;
      }
      const Box pin_box = toBox(pin.getBBox());
      LayerMetal metal = layerMetal(net, layer_num);
      for (const auto& group : metalOnPinComponents(metal.pieces, pin_box)) {
        components++;
        IoPinPatchRequest request;
        for (const size_t k : group) {
          request.group.push_back(metal.pieces[k]);
        }
        if (unionArea(request.group) >= area->getMinArea()) {
          continue;
        }
        below++;
        request.horizontal = layer->isHorizontal();
        request.half_width = layer->getWidth() / 2;
        request.min_area = area->getMinArea();
        request.grid = grid;
        request.die = die;
        // Try local anchors in survivor order; a rejected extension at one
        // anchor need not prevent a legal extension at another.
        for (const auto& anchor : metal.anchors) {
          if (!pin_box.containsClosed(anchor.x(), anchor.y())
              || !std::ranges::any_of(group, [&](const size_t k) {
                   return metal.pieces[k].containsClosed(anchor.x(), anchor.y());
                 })) {
            continue;
          }
          request.anchor = {anchor.x(), anchor.y()};
          request.legal = [&](const Box& candidate) {
            for (const Box& other : metal.pieces) {
              const GeoCoord prl = std::max<GeoCoord>(
                  0,
                  request.horizontal
                      ? std::min(candidate.xhi, other.xhi)
                            - std::max(candidate.xlo, other.xlo)
                      : std::min(candidate.yhi, other.yhi)
                            - std::max(candidate.ylo, other.ylo));
              const GeoCoord spacing = layer->getMinSpacingValue(
                  std::min(candidate.width(), candidate.height()),
                  std::min(other.width(), other.height()),
                  prl,
                  false);
              if (patchLeavesNarrowGap(candidate, other, spacing)) {
                return false;
              }
            }
            frPatchWire patch;
            patch.setLayerNum(layer_num);
            patch.setOrigin(anchor);
            odb::Rect offset(candidate.xlo, candidate.ylo, candidate.xhi, candidate.yhi);
            offset.moveDelta(-anchor.x(), -anchor.y());
            patch.setOffsetBox(offset);
            return patchClean(design, logger, router_cfg, net, patch);
          };
          const auto candidate = ioPinPatch(request);
          if (!candidate) {
            continue;
          }
          auto patch = std::make_unique<frPatchWire>();
          patch->setLayerNum(layer_num);
          patch->setOrigin(anchor);
          odb::Rect offset(candidate->xlo, candidate->ylo, candidate->xhi, candidate->yhi);
          offset.moveDelta(-anchor.x(), -anchor.y());
          patch->setOffsetBox(offset);
          frPatchWire* object = patch.get();
          net->addPatchWire(std::move(patch));
          design->getRegionQuery()->addDRObj(object);
          net->setModified(true);
          metal.pieces.push_back(*candidate);
          patched++;
          break;
        }
      }
    }
  }
  logger->info(utl::DRT,
                9263,
                "BonnRoute DR IO pins after FlexDR: {} routed components on "
                "top-level pin shapes, {} below the layer's minimum area "
                "without the pin, {} patched to it.",
                components,
                below,
                patched);
}

}  // namespace drt::br
