// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "br/BrSpeculativeRipup.h"

#include <omp.h>

#include <algorithm>
#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

#include "br/BonnRouteDR.h"

namespace drt {

void BonnRouteDR::searchRipupBatch(const std::vector<int>& queue,
                                   const size_t first,
                                   const int pass,
                                   RipupBatch& batch)
{
  batch.first = first;
  batch.previous_writes = {};
  batch.candidates.clear();
  const size_t end = std::min(queue.size(), first + params_.speculative_ripup);
  batch.candidates.resize(end - first);
  // Only these copies are mutable during presearch. The committed routing,
  // original nets and history remain fixed until all searches have ended.
  for (size_t k = first; k < end; k++) {
    const int id = queue[k];
    if (!nets_[id].routed) {
      batch.candidates[k - first].net = nets_[id];
      batch.candidates[k - first].net_version = batch.net_versions[id];
    }
  }
  const double begin = omp_get_wtime();
  const int threads = static_cast<int>(workers_.size());
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
  for (int k = 0; k < static_cast<int>(batch.candidates.size()); k++) {
    RipupCandidate& candidate = batch.candidates[k];
    if (!candidate.net) {
      continue;
    }
    Worker& local = workers_[omp_get_thread_num()];
    const br::AccessFootprintScope scope(&local.footprint);
    local.footprint.begin(candidate.net->id,
                          {0, 0, tiles_x_ - 1, tiles_y_ - 1});
    local.max_ripup_level
        = pass == 0
              ? static_cast<br::RipupLevel>(params_.first_pass_ripup_level)
              : br::kMaxRemovableLevel;
    local.defer_history = true;
    local.history_bumps.clear();
    const int64_t settled_start = local.settled;
    const double search_start = omp_get_wtime();
    candidate.found = connect(*candidate.net, true, local, pass);
    candidate.seconds = omp_get_wtime() - search_start;
    candidate.settled = local.settled - settled_start;
    candidate.reads = local.footprint.accessed();
    candidate.reads.add(local.footprint.wideReads());
    candidate.history = std::move(local.history_bumps);
    local.history_bumps.clear();
    local.defer_history = false;
  }
  batch.search_seconds += omp_get_wtime() - begin;
  for (const RipupCandidate& candidate : batch.candidates) {
    batch.searched += candidate.net.has_value();
    batch.searched_settled += candidate.settled;
  }
}

BonnRouteDR::RipupCandidate* BonnRouteDR::reuseRipupCandidate(
    const size_t index,
    const int id,
    RipupBatch& batch)
{
  RipupCandidate& candidate = batch.candidates[index - batch.first];
  if (!candidate.net) {
    return nullptr;
  }
  if (br::ripupSearchReusable(
          candidate.reads,
          batch.previous_writes,
          candidate.net_version == batch.net_versions[id])) {
    batch.reused++;
    batch.reused_settled += candidate.settled;
    return &candidate;
  }
  batch.invalidated++;
  return nullptr;
}

void BonnRouteDR::completeSerialRipup(const int id,
                                      const bool routed,
                                      const std::vector<SavedRouting>& journal,
                                      std::vector<int>& unrouted)
{
  if (routed) {
    return;
  }
  // Preserve the root failure before rollback restores the original routing.
  const uint8_t step = nets_[id].left_step;
  const int term = nets_[id].left_term;
  for (auto it = journal.rbegin(); it != journal.rend(); ++it) {
    restoreRouting(*it);
  }
  serial_ripup_.rollbacks++;
  unrouted.push_back(id);
  nets_[id].ripup_outcome = serial_ripup_.outcome_depth > 0
                                ? uint8_t{kRipupVictimFailed}
                                : serial_ripup_.outcome;
  nets_[id].left_step = step;
  nets_[id].left_term = term;
}

void BonnRouteDR::routeSerialRipupPass(std::vector<int>& queue,
                                       const int pass,
                                       const int64_t settled_budget,
                                       int64_t& serial_settled,
                                       RipupBatch* batch)
{
  Worker& worker = workers_.front();
  std::vector<int> unrouted;
  std::vector<SavedRouting> journal;
  if (batch != nullptr) {
    batch->candidates.clear();
    batch->first = 0;
  }
  for (size_t index = 0; index < queue.size(); index++) {
    const int id = queue[index];
    if (nets_[id].routed) {
      continue;
    }
    if (settled_budget > 0 && serial_settled >= settled_budget) {
      unrouted.push_back(id);
      ripup_budget_skipped_++;
      if (nets_[id].ripup_outcome == kRipupNotStarted) {
        nets_[id].ripup_outcome = kRipupBudget;
      }
      continue;
    }
    RipupCandidate* cached = nullptr;
    std::optional<br::AccessFootprint> writes;
    if (batch != nullptr) {
      if (index >= batch->first + batch->candidates.size()) {
        searchRipupBatch(queue, index, pass, *batch);
      }
      cached = reuseRipupCandidate(index, id, *batch);
      writes.emplace();
      writes->begin(id, {0, 0, tiles_x_ - 1, tiles_y_ - 1});
    }
    const br::AccessFootprintScope scope(
        writes ? &*writes : br::AccessFootprint::current());
    journal.clear();
    serial_ripup_.outcome = kRipupNotStarted;
    serial_ripup_.outcome_depth = 0;
    const int64_t sequence_start = worker.settled;
    const bool routed = ripupSequence(
        id, 0, pass, journal, false, worker, serial_ripup_, nullptr, cached);
    // Discarded searches consume resources, but never the serial budget.
    serial_settled += worker.settled - sequence_start
                      + (cached != nullptr ? cached->settled : 0);
    completeSerialRipup(id, routed, journal, unrouted);
    if (batch != nullptr) {
      // Even a rolled-back sequence changes metadata such as rip-up count.
      for (const SavedRouting& saved : journal) {
        batch->net_versions[saved.net]++;
      }
      for (int step = 0; step < br::AccessFootprint::kMaxSteps; step++) {
        batch->previous_writes.add(writes->writes(step));
      }
    }
  }
  queue.swap(unrouted);
}

}  // namespace drt
