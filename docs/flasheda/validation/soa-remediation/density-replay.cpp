#include <algorithm>
#include <chrono>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <vector>
#include "snapshot.h"
#include "legacyDensity.h"
#include "densitySoa.h"
using Clock = std::chrono::steady_clock;
int main(int argc, char** argv) {
  if (argc != 2) return 2;
  const auto s = flasheda::density::read(argv[1]);
  const int nx=s.nx, ny=s.ny;
  const size_t count=nx*ny;
  std::vector<int> lx(count),ly(count),ux(count),uy(count),xe(nx+1),ye(ny+1);
  std::vector<float> fx(count),fy(count),target(count);
  std::vector<int64_t> expected_inst(count),expected_fill(count),inst(count),fill(count);
  for (int y=0;y<ny;++y) for (int x=0;x<nx;++x) {
    const auto& b=s.bins[y*nx+x];const size_t i=x*ny+y;
    lx[i]=b.lx;ly[i]=b.ly;ux[i]=b.ux;uy[i]=b.uy;fx[i]=b.ex;fy[i]=b.ey;
    xe[x]=b.lx;xe[x+1]=b.ux;ye[y]=b.ly;ye[y+1]=b.uy;target[i]=b.target;
    expected_inst[i]=b.inst;expected_fill[i]=b.filler;
  }
  gpl::DensityGridView candidate{xe[0],ye[0],nx,ny,s.sx,s.sy,xe,ye,fx,fy};
  legacy_gpl::DensityGridView previous{xe[0],ye[0],nx,ny,s.sx,s.sy,lx,ly,ux,uy,fx,fy};
  // Resident cell columns are prepared once, outside every timed sample.
  std::vector<gpl::DensityRect> rects;std::vector<legacy_gpl::DensityRect> old_rects;
  for (const auto& c:s.cells) {rects.push_back({c.lx,c.ly,c.ux,c.uy,c.scale});old_rects.push_back({c.lx,c.ly,c.ux,c.uy,c.scale});}
  std::vector<std::pair<float,float>> gradients(s.cells.size());
  auto run=[&](auto& grid,auto& cells,int repeats) {
    double scatter=0,gather=0;
    for(int r=0;r<repeats;++r) {
      auto start=Clock::now();
      std::fill(inst.begin(),inst.end(),0);std::fill(fill.begin(),fill.end(),0);
      for(size_t i=0;i<cells.size();++i) {
        const int kind=s.cells[i].kind;if(!kind)continue;
        auto& out=kind==3?fill:inst;
        grid.forEachOverlap(cells[i],[&](size_t b,float area){out[b]+=static_cast<int64_t>(kind==2?area*target[b]:area);});
      }
      scatter+=std::chrono::duration<double,std::milli>(Clock::now()-start).count();
      if(inst!=expected_inst || fill!=expected_fill)throw std::runtime_error("captured scatter differs");
      start=Clock::now();
      for(size_t i=0;i<cells.size();++i)gradients[i]=grid.gather(cells[i]);
      gather+=std::chrono::duration<double,std::milli>(Clock::now()-start).count();
      for(size_t i=0;i<cells.size();++i)if(gradients[i].first!=s.cells[i].gx || gradients[i].second!=s.cells[i].gy)throw std::runtime_error("captured gather differs");
    }
    return std::pair{scatter/repeats,gather/repeats};
  };
  run(previous,old_rects,1);run(candidate,rects,1);
  for(int sample=0;sample<11;++sample) {
    for(int order=0;order<2;++order) {
      const bool old=(order+sample)%2==0;
      auto [scatter,gather]=old?run(previous,old_rects,20):run(candidate,rects,20);
      std::cout << "{\"snapshot\":\"" << argv[1] << "\",\"variant\":\"" << (old?"previous":"separable") << "\",\"sample\":" << sample << ",\"scatter_ms\":" << scatter << ",\"gather_ms\":" << gather << "}\n";
    }
  }
}
