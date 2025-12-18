#include "TASImage.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TError.h"
#include "TStopwatch.h"
#include <vector>
#include <algorithm>
#include <cstdint>
#include <cmath>
#include <iostream>

/*
This script takes a source image and a target image to create a map from the source pixels to reconstruct the target image
To do this, all of the pixels from each image are ordered by RGB in base 256 and an initial map of ordered source to ordered target
is chosen. Energy is the difference in RGB from the map to the target (defined like physical distance). The optimization process is
swapping two pixels randomly (within a certain RADIUS to avoid poor swaps) and seeing how the energy changes. The metropolis algorithm 
is used for accepting changes and this process is done as simulated annealing with exponential decay of temperature.
*/

//setting up the source and target images as well as the output
//so just flip image1 and imag2 to get the other output
static const char* SRC_IMAGE = "./image1.png";
static const char* TGT_IMAGE = "./image2.png";
static const char* OUT_IMAGE = "./out1to2.png";

//parameters
static const int N_STEPS = 40000000;
static const int RADIUS = 8000; //we order by color then do swaps. restrict the radius so we don't try obvious bad swaps
static const double T_I = 2.0; //initial temp 
static const double T_F = 1e-3; //final temp


// Pixel helpers (ARGB = 0xAARRGGBB)
inline int R(UInt_t p) { return (p >> 16) & 0xFF; }
inline int G(UInt_t p) { return (p >> 8) & 0xFF; }
inline int B(UInt_t p) { return (p) & 0xFF; }

//treat RGB like xyz so we can do our normal distance formula. The further two colors are, the higher the energy (worse match)
inline double pixEnergy(UInt_t src, UInt_t tgt) {
  const double dr = R(src) - R(tgt);
  const double dg = G(src) - G(tgt);
  const double db = B(src) - B(tgt);
  return dr*dr + dg*dg + db*db; 
}

//colors in RGB go from 0-255, so we can write all the info of one pixel's color as RGB as an int
//it's like storing bits where the 1's place is B, 10's place is G, and 100's place is R (base 256)
inline uint64_t colorKey(UInt_t p) {
  const uint64_t r = (uint64_t)R(p);
  const uint64_t g = (uint64_t)G(p);
  const uint64_t b = (uint64_t)B(p);
  return r * 256 * 256 + g * 256 + b; //weighting to that place splitting thing
}

//we don't want to forget what initial position is associated with what pixel
//the pixels are gonna be sorted by color, so we need to define how one key is less than another
struct KeyIndex {
  uint64_t key;
  int idx;

  bool operator<(const KeyIndex& other) const {
    return key < other.key; 
  }
};

//temp decrease for annealing
double temperature(int step) {
  double fraction = double(step) / double(N_STEPS - 1);
  double T = T_I * std::pow(T_F / T_I, fraction);
  return T;
}

void pixswap() {

  TStopwatch timer;
  timer.Start();

  //load the images. we set them earlier so we don't gotta change this for swapping the images
  TASImage srcImg(SRC_IMAGE);
  TASImage tgtImg(TGT_IMAGE);

  int w = srcImg.GetWidth(); //1920, could've hardcoded it 
  int h = srcImg.GetHeight(); //1080
  int N = w * h; //number of pixels so 1920*1080=2073600
  
  //RGB array. This includes alpha which is transparency but we don't care about that
  UInt_t* srcPix = srcImg.GetArgbArray();
  UInt_t* tgtPix = tgtImg.GetArgbArray();

  //makes vectors (so arrays lol) of the key index which is the RGB and the initial position 
  std::vector<KeyIndex> srcK(N);
  std::vector<KeyIndex> tgtK(N);
  for (int i = 0; i < N; ++i) {
    //setting the RGB and initial position. Since we list all pixels at once, their ith index tells the exact position
    srcK[i] = { colorKey(srcPix[i]), (int)i }; 
    tgtK[i] = { colorKey(tgtPix[i]), (int)i };
  }

  //organize the pixels by RGB
  std::sort(srcK.begin(), srcK.end());
  std::sort(tgtK.begin(), tgtK.end());

  //list the original positions now that they're sorted by RGB (so [0,1,2,3] may have become [2,1,3,0] when sorted)
  std::vector<int> srcOrder(N);
  std::vector<int> tgtOrder(N);
  for (int i = 0; i < N; ++i) {
    srcOrder[i] = srcK[i].idx;
    tgtOrder[i] = tgtK[i].idx;
  }

  //doing a random start is a terrible idea
  //so we start with the map just being the best line up by RGB sorting
  std::vector<int> mapping(N);
  for (int i = 0; i < N; ++i)
    mapping[tgtOrder[i]] = srcOrder[i];

  //initial energy
  double E = 0.0;
  for (int i = 0; i < N; ++i)
    E += pixEnergy(srcPix[mapping[i]], tgtPix[i]);

  TRandom3 rng(12345);

  //Simulated Annealing
  //swap two pixels' source pixel colors and check how that changes the energy
  //don't recompute all energy every time since thats a time sink
  for (int step = 0; step < N_STEPS; ++step) {
    const double T = temperature(step);

    int k = rng.Integer((UInt_t)N);
    int dk = rng.Integer(2*RADIUS + 1) - RADIUS;
    if (dk == 0) dk = 1;

    int k2 = k + dk;
    if (k2 < 0) k2 = 0;
    if (k2 >= (int)N) k2 = (int)N - 1;
    if (k2 == k) continue;

    int i = tgtOrder[k];
    int j = tgtOrder[k2];

    int mi = mapping[i];
    int mj = mapping[j];

    //dE = E_f - E_i -> (Energy for map(j) with i + map(i) with j) - (Energy for map(i) with i + energy for map(j) with j)
    double dE = pixEnergy(srcPix[mj], tgtPix[i]) + pixEnergy(srcPix[mi], tgtPix[j]) - 
        pixEnergy(srcPix[mi], tgtPix[i]) - pixEnergy(srcPix[mj], tgtPix[j]);

    //metropolis algorithm
    if (dE <= 0.0 || rng.Rndm() < std::exp(-dE / T)) {
      mapping[i] = mj;
      mapping[j] = mi;
      E += dE;
    }
  }


  //take our "optimized" map and make the image
  TASImage outImg(srcImg);
  UInt_t* outPix = outImg.GetArgbArray();

  for (int i = 0; i < N; ++i)
    outPix[i] = srcPix[mapping[i]];

  outImg.WriteImage(OUT_IMAGE);

  timer.Stop();
  std::cout << "Total runtime = " << timer.RealTime() << " s\n";  //Total runtime = 15.8127 s from one of the runs
}
