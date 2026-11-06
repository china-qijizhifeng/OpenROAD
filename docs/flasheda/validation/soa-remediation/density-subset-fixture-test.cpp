// Validation-only adapter over the existing integrated placement fixture.
#include <bit>
#include <cstdint>
#include <stdexcept>

#include "src/gpl/src/backendContext.h"
#include "src/gpl/src/densityGradientBackend.h"
#include "src/gpl/test/congestion_buffers_timing_test.cpp"

namespace gpl {
TEST_F(TimingRepairIntegrationTest, DensityBackendFullSubsetAndReleaseSizes)
{
  makePlacer(false, false);
  for (const auto& nb : nb_vec_) {
    ASSERT_FALSE(nb->getGCells().empty());
    auto backend = makeDensityGradientBackend(BackendContext{.nb = nb.get()});
    const auto& cells = nb->getGCells();
    PointColumns all;
    all.resize(cells.size());
    backend->getCellGradients(cells, all);
    for (size_t i = 0; i < cells.size(); ++i) {
      const auto scalar = nb->getDensityGradient(cells[i]);
      EXPECT_EQ(std::bit_cast<uint32_t>(all.x()[i]),
                std::bit_cast<uint32_t>(scalar.x));
      EXPECT_EQ(std::bit_cast<uint32_t>(all.y()[i]),
                std::bit_cast<uint32_t>(scalar.y));
    }
    std::vector<GCellHandle> subset;
    for (size_t i = cells.size(); i > 0; i -= std::min<size_t>(i, 3)) {
      subset.push_back(cells[i - 1]);
    }
    subset.push_back(subset.front());
    PointColumns out;
    out.resize(subset.size());
    backend->getCellGradients(subset, out);
    for (size_t i = 0; i < subset.size(); ++i) {
      const auto scalar = nb->getDensityGradient(subset[i]);
      EXPECT_EQ(std::bit_cast<uint32_t>(out.x()[i]),
                std::bit_cast<uint32_t>(scalar.x));
      EXPECT_EQ(std::bit_cast<uint32_t>(out.y()[i]),
                std::bit_cast<uint32_t>(scalar.y));
    }
    std::vector<GCellHandle> empty;
    PointColumns empty_out;
    EXPECT_NO_THROW(backend->getCellGradients(empty, empty_out));
    PointColumns wrong;
    wrong.resize(subset.size() + 1);
    EXPECT_THROW(backend->getCellGradients(subset, wrong), std::invalid_argument);
    EXPECT_THROW(backend->getCellGradients(cells, empty_out),
                 std::invalid_argument);
  }
}
}  // namespace gpl
