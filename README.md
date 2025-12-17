# simpix

C++ starter code
* simpix_start.cpp
use make to build this example

Usage: simapix_start image1 image2 <output=out.png>

Python starter code
* simpix_start.py

Usage: simapix_start image1 image2 <output=out.png>


pixswap.cpp - this is a ROOT macro so execute as .x pixswap.cpp in ROOT. I left comments throughout to hopefully make it clear how it runs

At the top of the script, you define your source and target images, so I got both images by just flipping image1 and image2 in these two lines. image1.png and image2.png are my original images, and out1to2.png and out2to1.png are the outputs from the script. Image size is not hardcoded so in theory, this script can be run on other size images. However, I had to do a smarter approach to test permutations, so you need to watch out for the variable RADIUS

Some smart decisions:
1. Don't initialize with a completely random map. The original map is just "the ith pixel by RGB ordering for the source image is the ith pixel for the target image". In theory this starting point should be pretty close to our goal (reds, greens, blues should be well aligned)
2. different configurationsa are tested by swapping two pixels: limit how far two pixels can be before they are swapped (RADIUS). Clearly if a section of the image is dominated by shades of red, trying out a random green pixel isn't smart, so our swaps should be two red-adjacent pixels instead. RADIUS has a hardcoded value, so if you try a smaller image you probably need to change this.

This code only takes around 15 seconds to run for my 1920x1080 images (Total runtime = 15.8127 s and Total runtime = 15.7577 s). What cuts down on time is not recomputing the entire energy after each swap. Instead the script just looks at how energy changes from the swap. For every pixel but the two swapped ones, the energy is unchanged so we don't want to compute those ones. We just want to know the difference in energy for those two pixels 