Author: Emmanuel Maggiori. March 2014.

Complimentary material for the literature review:
"Perceptual grouping by tensor voting: a comparative survey of recent approaches". E Maggiori, HL Manterola, M del fresno. To be published in IET Computer Vision.

Implementation of Steerable Tensor Voting as published in:
"An efficient method for tensor voting using steerable filters", Franken et. al. ECCV 2006.

Usage (example):

[s,o]=encode('Lena.png'); %To encode the tensorized gradient of an image
[saliency,ballness,orientation]=vote(s,o,5); %Tensor Voting with sigma=5



