This is the toolbox V1.0 to measure the VTs with human subjects. Author: Feng Liu, University of Arizona, Tucson, AZ, USA, liuf@email.arizona.edu

To do the Visual Threshold (VT) measurements, please run the "runme.m" manuscript in Matlab.

For each of the luminance, Cb and Cr components, measurements of VTs at lower JND levels has to be totaly finished before VTs at higher JND levels are measured.

VTs for HH1-5, HL1-5, LH1-5 and LL5 are needed at each JND level for the encoder. The VTs of HL and LH subbands at a JND level are very close to each other. So the final results for these two subbands are pooled by taking the minimum. 

For a non-LL subband of the luminance component, we model the VT versus the subband coefficient variance linearly. To find the linearity parameters (slope and intersection), VT measurements for 4 subband coefficient variances are needed for each subband. For non-LL subbands of luminance component, HL and LH subbands at a given decomposition level share the same coefficient variances. So at a given decomposition level of luminance component, please finish the HL subband measurements before proceeding to the LH subband. For LL5 subband of luminance component, only 1 coefficient variance value is needed.

For a subband of chrominance component, the VT do not vary a lot as the subband coefficient variance changes. So a single typical variance from [1] is used for each subband.

The measurement procedure is based on the QUEST psychometric method [2]. Part of the files in the toolbox (those beginning with "Quest") are refered to [3]. Please also referred to the license of the QUEST software in "./toolbox32/".

3AFC trials with the question "Which stimulus is different from the other two?" are taken for each subband at a given subband coefficient variance and JND level. 10-second observation time limit is enabled for each trial. The subject has to answer the question by pressing "1", "2" or "3" at the end of each trial.

For more details of the measurement experiment set-up and procedure, as well as how to use the VTs in the encoder, please refer to [4].

The measured results for luminance component will be stored in "./toolbox/lumthresh/", and measured results for chrominance components will be stored in "./toolbox/chromthresh/". The data to be read by the encoder will be stored in "./toolbox/data/".

[1] H. Oh, A. Bilgin and M. W. Marcellin, "Visually Lossless Encoding for JPEG2000," in IEEE Transactions on Image Processing, vol. 22, no. 1, pp. 189-201, January, 2013.

[2] A. B. Watson and D. G. Pelli, "QUEST: A Bayesian adaptive psychometirc method," Perception and psychophysics, vol. 33, no. 2, February, 1983.

[3] "Psychtoolbox-3 software [online].", http://psychtoolbox.org/.

[4] Feng Liu, Yuzhang Lin, Eze L. Ahanonu, Michael W. Marcellin, Amit Ashok, Elizabeth A. Krupinski and Ali Bilgin, "Visibility Thresholds for Visually Lossy JPEG2000," Proc. SPIE 9971, Applications of Digital Image Processing XXXIX, 99711P, September, 2016. 