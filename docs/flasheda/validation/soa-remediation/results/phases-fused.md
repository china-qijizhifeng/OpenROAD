| phase | fused-1 | fused-8 |
|---|---:|---:|
| timing_feedback | 17.81 | 17.86 |
| routability_feedback | 0.81 | 0.76 |
| initial_placement | 3.79 | 1.51 |
| density_scatter | 9.78 | 2.49 |
| fft | 5.02 | 1.34 |
| density_field_other | 0.55 | 0.09 |
| density_gather | 10.71 | 1.60 |
| wirelength_force | 2.49 | 0.51 |
| wirelength_gather | 3.38 | 1.96 |
| gradient_combine | 1.23 | 0.73 |
| hpwl | 0.90 | 0.10 |
| coordinates | 4.76 | 0.73 |
| step_length | 1.61 | 0.93 |
| update_iteration | 0.72 | 0.17 |
| nesterov_other | 0.04 | 0.04 |
| gpl_other | 0.33 | 0.33 |
| **total** | **63.91** | **31.16** |

### timing_feedback
- fused-1: sta::Table::findValueOrder2 const 0.80s; sta::exp2 0.50s; sta::findValueIndex 0.43s; odb::dbITerm::getInst const 0.38s; sta::DmpAlg::Vl 0.38s; sta::findRoot 0.36s; std::_Rb_tree::_M_get_insert_unique_pos 0.34s; sta::dbNetwork::staToDb const 0.29s
- fused-8: sta::Table::findValueOrder2 const 0.81s; sta::exp2 0.47s; sta::findRoot 0.42s; sta::DmpAlg::Vl 0.39s; std::_Rb_tree::_M_get_insert_unique_pos 0.38s; sta::findValueIndex 0.37s; odb::dbITerm::getInst const 0.34s; std::vector const& est::EstimateParasitics::resolveWireRC const 0.30s

### routability_feedback
- fused-1: grt::Rudy::processIntersectionSignalNet 0.10s; gpl::BinGrid::scatterDensityAreaInPlace 0.06s; [libc.so.6] <- boost::multi_array::resize 0.04s; grt::Rudy::Tile::addDemand 0.04s; gpl::NesterovBase::updateDensitySize 0.04s; odb::dbITerm::getMTerm const 0.04s; gpl::NesterovBase::getFillerGCell 0.02s; odb::dbITerm::getAvgXY const 0.02s
- fused-8: grt::Rudy::processIntersectionSignalNet 0.11s; [libc.so.6] <- boost::multi_array::resize 0.04s; odb::dbITerm::getMTerm const 0.04s; grt::Rudy::Tile::addDemand 0.03s; odb::dbTransform::apply const 0.03s; grt::FastRouteCore::getCapacityReductionData 0.02s; gpl::BinGrid::scatterDensityAreaInPlace 0.02s; odb::dbNet::getName[abi:cxx11] const 0.02s

### initial_placement
- fused-1: bool Eigen::internal::bicgstab 3.04s; gpl::InitialPlace::createSparseMatrix 0.11s; void Eigen::internal::set_from_triplets 0.08s; odb::dbITerm::getInst const 0.05s; [libc.so.6] <- bool Eigen::internal::bicgstab 0.05s; Eigen::SparseMatrix& Eigen::SparseMatrix::operator= 0.04s; gpl::Pin::cx const 0.04s; odb::dbInst::getBBox 0.03s
- fused-8: Eigen::internal::sparse_time_dense_product_impl::run 0.47s; bool Eigen::internal::bicgstab 0.17s; [libc.so.6] <- bool Eigen::internal::bicgstab 0.12s; gpl::InitialPlace::createSparseMatrix 0.12s; void Eigen::internal::set_from_triplets 0.08s; [libgomp.so.1.0.0] <- bool Eigen::internal::bicgstab 0.07s; Eigen::SparseMatrix& Eigen::SparseMatrix::operator= 0.05s; odb::dbITerm::getInst const 0.04s

### density_scatter
- fused-1: gpl::BinGrid::scatterDensityAreaInPlace 7.99s; gpl::BinGrid::updateBinsGCellDensityArea 1.71s; [libc.so.6] <- gpl::BinGrid::updateBinsGCellDensityArea 0.08s; odb::dbDatabase::getChip 0.00s; odb::dbChip::getBlock 0.00s
- fused-8: gpl::BinGrid::scatterDensityAreaInPlace 1.31s; gpl::BinGrid::updateBinsGCellDensityArea 0.88s; [libc.so.6] <- gpl::BinGrid::updateBinsGCellDensityArea 0.10s; [libc.so.6] <- gpl::BinGrid::scatterDensityAreaInPlace 0.09s; [libgomp.so.1.0.0] <- gpl::BinGrid::updateBinsGCellDensityArea 0.06s; [libgomp.so.1.0.0] <- gpl::BinGrid::scatterDensityAreaInPlace 0.05s; gpl::BinGrid::densityView const 0.00s; odb::dbDatabase::getChip 0.00s

### fft
- fused-1: gpl::cftf1st 0.59s; gpl::::CpuFftBackend::solve 0.52s; gpl::dctsub 0.47s; gpl::ddxt2d_sub 0.45s; gpl::rftfsub 0.45s; gpl::cftmdl1 0.45s; gpl::cftf081 0.40s; gpl::ddct 0.29s
- fused-8: gpl::::CpuFftBackend::solve 0.38s; gpl::::CpuFftBackend::transform 0.18s; [libgomp.so.1.0.0] <- gpl::::CpuFftBackend::transform 0.18s; gpl::cftf1st 0.07s; [libc.so.6] <- gpl::FFT::solveColumns 0.07s; gpl::cftmdl1 0.06s; gpl::cftf081 0.06s; gpl::dctsub 0.06s

### density_field_other
- fused-1: gpl::NesterovBase::updateDensityFieldBin 0.55s; [libgomp.so.1.0.0] <- gpl::NesterovBase::updateDensityFieldBin 0.00s; [libc.so.6] <- gpl::NesterovBase::updateDensityFieldBin 0.00s
- fused-8: gpl::NesterovBase::updateDensityFieldBin 0.08s; [libgomp.so.1.0.0] <- gpl::NesterovBase::updateDensityFieldBin 0.01s

### density_gather
- fused-1: void gpl::DensityGridView::forEachOverlap const 7.86s; void gpl::NesterovBase::combineGradientsColumns 2.85s; [libc.so.6] <- void gpl::NesterovBase::combineGradientsColumns 0.00s; [libgomp.so.1.0.0] <- void gpl::NesterovBase::combineGradientsColumns 0.00s
- fused-8: void gpl::DensityGridView::forEachOverlap const 1.05s; void gpl::NesterovBase::combineGradientsColumns 0.41s; [libgomp.so.1.0.0] <- void gpl::NesterovBase::combineGradientsColumns 0.14s

### wirelength_force
- fused-1: gpl::updateNetForceColumns 2.37s; [libc.so.6] <- gpl::NesterovBaseCommon::updateWireLengthForceWA_soa 0.07s; gpl::NesterovBaseCommon::updateWireLengthForceWA_soa 0.05s; [libgomp.so.1.0.0] <- gpl::NesterovBaseCommon::updateWireLengthForceWA_soa 0.00s
- fused-8: gpl::updateNetForceColumns 0.39s; [libc.so.6] <- gpl::NesterovBaseCommon::updateWireLengthForceWA_soa 0.08s; [libgomp.so.1.0.0] <- gpl::NesterovBaseCommon::updateWireLengthForceWA_soa 0.03s; gpl::NesterovBaseCommon::updateWireLengthForceWA_soa 0.01s

### wirelength_gather
- fused-1: gpl::NesterovBaseCommon::getWireLengthGradientWA_soa const 2.95s; gpl::GCell::isFiller const 0.43s
- fused-8: gpl::NesterovBaseCommon::getWireLengthGradientWA_soa const 1.92s; gpl::GCell::isFiller const 0.04s

### gradient_combine
- fused-1: gpl::NesterovBase::updateGradients 0.67s; gpl::NesterovBaseCommon::getWireLengthGradientsWA_cpu 0.27s; gpl::NesterovBase::getFillerGCell 0.25s; gpl::NesterovBaseCommon::getGCell 0.04s; [libgomp.so.1.0.0] <- gpl::NesterovBaseCommon::getWireLengthGradientsWA_cpu 0.00s
- fused-8: gpl::NesterovBase::updateGradients 0.65s; gpl::NesterovBaseCommon::getWireLengthGradientsWA_cpu 0.05s; gpl::NesterovBaseCommon::getGCell 0.03s; [libgomp.so.1.0.0] <- gpl::NesterovBaseCommon::getWireLengthGradientsWA_cpu 0.01s; gpl::NesterovBase::getFillerGCell 0.00s

### hpwl
- fused-1: gpl::NesterovBaseCommon::getHpwlColumns 0.90s; [libgomp.so.1.0.0] <- gpl::NesterovBaseCommon::getHpwlColumns 0.00s
- fused-8: gpl::NesterovBaseCommon::getHpwlColumns 0.09s; [libgomp.so.1.0.0] <- gpl::NesterovBaseCommon::getHpwlColumns 0.01s; gpl::NesterovBaseCommon::getHpwl 0.00s

### coordinates
- fused-1: gpl::NesterovBase::updateNextCoordinatesColumns 2.83s; gpl::NesterovBase::updateGCellDensityCenterLocation 0.94s; gpl::NesterovBaseCommon::updateDensityPinLocations 0.52s; gpl::BinGrid::ux const 0.16s; gpl::BinGrid::lx const 0.11s; gpl::BinGrid::uy const 0.10s; gpl::BinGrid::ly const 0.10s; [libgomp.so.1.0.0] <- gpl::NesterovBase::updateNextCoordinatesColumns 0.00s
- fused-8: gpl::NesterovBase::updateNextCoordinatesColumns 0.36s; gpl::NesterovBase::updateGCellDensityCenterLocation 0.14s; gpl::NesterovBaseCommon::updateDensityPinLocations 0.11s; [libgomp.so.1.0.0] <- gpl::NesterovBase::updateNextCoordinatesColumns 0.03s; [libgomp.so.1.0.0] <- gpl::NesterovBase::updateGCellDensityCenterLocation 0.03s; gpl::BinGrid::ux const 0.02s; gpl::BinGrid::lx const 0.02s; [libgomp.so.1.0.0] <- gpl::NesterovBaseCommon::updateDensityPinLocations 0.02s

### step_length
- fused-1: gpl::getDistance 1.60s; [libgomp.so.1.0.0] <- gpl::getDistance 0.01s; std::vector::_M_fill_insert 0.00s; [libc.so.6] <- gpl::getDistance 0.00s
- fused-8: gpl::getDistance 0.91s; [libgomp.so.1.0.0] <- gpl::getDistance 0.02s; std::vector::_M_fill_insert 0.00s

### update_iteration
- fused-1: gpl::GCell::isInstance const 0.31s; gpl::NesterovBase::getFillerGCell 0.16s; gpl::NesterovBase::updateNextIter 0.14s; gpl::NesterovBaseCommon::getGCell 0.07s; gpl::GCell::isLocked const 0.03s; [libc.so.6] <- spdlog::sinks::ansicolor_sink::log 0.00s; [libm.so.6] <- gpl::NesterovBase::getPhiCoef const 0.00s; char const* fmt::v11::detail::parse_format_specs 0.00s
- fused-8: gpl::NesterovBaseCommon::getGCell 0.06s; gpl::GCell::isInstance const 0.04s; gpl::GCell::isLocked const 0.03s; gpl::NesterovBase::updateNextIter 0.03s; [libgomp.so.1.0.0] <- gpl::NesterovBase::updateNextIter 0.00s; gpl::NesterovBase::getFillerGCell 0.00s; [libc.so.6] <- spdlog::sinks::ansicolor_sink::log 0.00s; [libm.so.6] <- gpl::NesterovBase::getPhiCoef const 0.00s

### nesterov_other
- fused-1: grt::Rudy::processIntersectionSignalNet 0.01s; odb::dbITerm::getMTerm const 0.01s; grt::Rudy::Tile::addDemand 0.00s; grt::FastRouteCore::getCapacityReductionData 0.00s; odb::dbITerm::getInst const 0.00s; odb::dbITerm::getAvgXY const 0.00s; [libc.so.6] <- std::vector::operator= 0.00s; odb::dbTransform::apply const 0.00s
- fused-8: grt::Rudy::processIntersectionSignalNet 0.01s; grt::Rudy::Tile::addDemand 0.01s; odb::dbInst::getTransform 0.00s; odb::dbTransform::apply const 0.00s; odb::dbITerm::getMTerm const 0.00s; odb::dbITerm::getAvgXY const 0.00s; [libc.so.6] <- std::vector::operator= 0.00s; odb::dbInst::getMaster const 0.00s

### gpl_other
- fused-1: cuddInitCache 0.03s; sta::Levelize::findTopologicalOrder 0.02s; gpl::SoaRow::SoaRow 0.01s; odb::dbITerm::getMTerm const 0.01s; gpl::NesterovBaseCommon::NesterovBaseCommon 0.01s; int& std::vector::emplace_back 0.01s; std::vector::_M_fill_insert 0.01s; odb::dbInstITermItr::next const 0.01s
- fused-8: cuddInitCache 0.04s; sta::Levelize::findTopologicalOrder 0.03s; gpl::SoaRow::SoaRow 0.01s; int& std::vector::emplace_back 0.01s; gpl::NesterovBaseCommon::NesterovBaseCommon 0.01s; [libc.so.6] <- cuddFreeTable 0.01s; sta::Graph::makeEdge 0.01s; [libc.so.6] <- cuddInitTable 0.01s
