// Temporary targeted connectivity diagnostic; removed after diagnosis.
#pragma once

#include <cstdlib>
#include <sstream>
#include <string>

#include "db/obj/frBTerm.h"
#include "db/obj/frInstTerm.h"
#include "db/obj/frNet.h"
#include "db/obj/frShape.h"
#include "db/obj/frVia.h"
#include "utl/Logger.h"

namespace drt {
inline bool connectivityProbeSelected(const frNet* net)
{
  const char* selected = std::getenv("DRT_CONNECTIVITY_PROBE_NET");
  return selected != nullptr && net->getName() == selected;
}

inline void connectivityProbe(const frNet* net,
                              const char* phase,
                              int iter,
                              utl::Logger* logger)
{
  if (!connectivityProbeSelected(net)) {
    return;
  }
  std::ostringstream out;
  out << "[DEBUG-c9300] net " << net->getName() << " phase " << phase
      << " iter " << iter << " initial " << net->hasInitialRouting() << '\n';
  for (const auto& obj : net->getShapes()) {
    out << *obj << '\n';
  }
  for (const auto& obj : net->getVias()) {
    out << *obj << " connected " << obj->isBottomConnected() << ' '
        << obj->isTopConnected() << '\n';
  }
  for (const auto& obj : net->getPatchWires()) {
    out << *obj << '\n';
  }
  for (const auto* term : net->getInstTerms()) {
    out << *term << '\n';
  }
  for (const auto* term : net->getBTerms()) {
    out << *term << '\n';
  }
  logger->report("{}", out.str());
}
}  // namespace drt
