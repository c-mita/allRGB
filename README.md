# All RGB Nearest Neighbour filling

A program to generate images containing all RGB values, incrementally
according to a nearest neighbour algorithm.

Primary implemneted using a KD-Tree, but a Vantage-Point (VP) Tree is also
implemented that can use the L-1, L-2 or L-infinity metrics.

## Build
```
$ zig build --release=fast
```

## Usage
```
$ ./zig-out/bin/main --help
...
$ ./zig-out/bin/main
# generates out.png
```

By default a full 24-bit 4096x4096 image is built using a KD-tree with pixels
initially generated in hue order.

For faster image generation the bit depth can be lowered:
```
$ ./zig-out/bin/main --depth 7 --out out.png
```
Will generate a 21 bit image (2048x1024 pixels) and is considerably faster.

Additionally a source image can be provided which will be used as the source
of pixels (including duplication). The output image will match the size.

## Implementation notes

The primary data structure is a KD-Tree with some alterations. Non-leaf nodes
do not correspond to any particular key (pixel value) and leaf nodes are buckets,
containing many (128) permitting a fast linear scan.

During image generation, a pixel is pulled from the source list and its
nearest-neigbour in the Tree is retrieved. Pixels in the tree correspond to
pixels in the output image that still have open neighbours (pixels that have
yet to be populated).
When a pixel has its last open neighbour removed it is removed from the tree.

KD-Trees do not naturally support key removal - the use of buckets makes this
more palatable but since splitting planes are determined by pixels that may
have all been removed and do not well fit the curret set under consideration,
the tree structure naturally degrades over the course of image generation.

To handle this we heuristically scrap the current tree and rebuild it from the
partially populated output image.

## Approximate matching

The `--approx` flag permits the program to deviate slightly from
nearest-neighbour.

For a given node in the tree, the side containing the query key is searched
for a match and if one is found this is simply returned. Only if no key is
found will the other side of the tree but examined.

Depending on the source and ordering, this can generate result in either minor
or major differences in the output.

## VP Trees

Optionally a Vantage-Point Tree can be used. These are more general than
KD or BSP trees since they require only a metric (instead of something like
an inner-product vector space) and offer flexibility in how distances are
calculated.

However they generally perform a little worse for this task when the true
nearest-neighbour algorithm is used.

With the `--approx` flag VP trees (particularly with the L1 metric) seem to
produce images closer to `--noapprox` than the KD tree does.


## Zig version

Targets and compiles with zig 0.16
