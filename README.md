# Supernova Classification by Host Galaxy with Machine Learning
### Center for Astrophysics | Harvard and Smithsonian 

_Formatting for publication in proccess_ <p>
This software classifies supernovae by type (Ia, Ibc, II, IIn, or Superluminous) using archival images of the supernova's host galaxy. \n
The main scripts are  
Host galaxy detection and property extraction: src/host_property_extractor.py  
Random forrest classifier: src/random_forrest_classifier.py  
Convolutional Neural Network (in progress): src/cnn_classifier.py  

## Description
Optical surveys are discovering supernovae at exponentially increasing rates. In 2023, LSST (the Large Synoptic Survey Telescope) will begin surveying, and it is predicted to discover 1 million supernovae per year. The current rate of discovery is about 10,000 per year, so this is a 100-fold increase, the largest increase in supernova detection in history. The traditional method of classifying supernovae using spectroscopy is expensive and slow, so only a fraction of these discovered supernovae will receive spectroscopic classification. Therefore, developing other classification methods is currently important and urgent. <p>
A great amount of work has been done in developing classification methods using photometry. However, little has been done using the host galaxy of a supernova, which also contains valuable information for predicting the type of the supernova. <p>
My overall goal was to classify supernova based on the host galaxy in which they occur. More specifically, my goals were to develop two machine learning algorithms accomplishing this. I trained both algorithms using about 400 spectroscopically classified supernovae from the Pan-STARRS medium-deep survey. My sample consisted of types Ia, Ibc, II, IIn, and superluminous supernovae, and my goal was to classify between those types. I will explain these two algorithms and describe my challenges and success in each. <p>
My first algorithm used preset galaxy properties, such as radius, magnitude (brightness), and ellipticity. It took images of the locations where the supernovae occurred, detected the host galaxy, and extracted properties from it. Correctly identifying the host galaxy proved a challenge I struggled with for most of the summer. About 80% of the time, there is an obvious host galaxy in the image, but the rest of the time there are complications. Because we see a two-dimensional projection of the sky, sometimes there are multiple galaxies that appear to be close to the supernova location, while actually some of them are actually far from the supernova because they are much closer to earth but just happen to be along a similar line of sight. I addressed this by looking up and take into account the redshift (which is correlated with distance from earth) of the objects in the image. Furthermore, sometimes the best host galaxy candidate would show up in one light wavelength but be too dim in another, so I cross-referenced across images taken in different wavelengths. Other times a close-by star, being very bright in the image, would get in the way of detection, so I used a catalog to identify and ignore the stars. In cases of star interference, I also set a higher threshold for the brightness necessary for my algorithm to pay attention to pixels, in order for it not to get confused by the way the star brightened a large section of the image. In the end, I was able to get my algorithm to detect a good host in about 90% of the images. The rest were mostly cases where there was clearly no detectable host galaxy in the image, plus a few very unusual cases mixed in. This success rate was sufficient for my purposes. <p>
I then fed these extracted properties into a random forest classifier, which is a machine learning algorithm consisting of a collection of decision trees that collectively classify an input. My main challenge in training the random forest classifier was the small number of samples I had that belonged to rarer classes of supernovae. I only had 9 samples of Superluminous supernovae, for example, as compared with over 200 of Type Ia supernovae. Firstly, the low number of the rarer classes makes it difficult for the machine to learn robust patterns that would be characteristic of that type. Furthermore, the imbalanced proportions of my sample make the machine tend to classify everything as the majority type, since by doing so it gets a high success rate. I addressed this by using SMOTE oversampling. This is a technique for generating synthetic data, which means a collection of artificial data points which resemble the original dataset. There are accuracy drawbacks to this solution, and it still did not perform nearly as well as a larger balanced input sample would. With it, though, I was still able to achieve about 45% average accuracy per class among 5 classes. For comparison, random classification would have yielded 20% accuracy. <p>
My second algorithm used a convolutional neural network. I created a basic form of this algorithm over the summer, and I plan to continue to improve and develop it during the semester. I trained it using the same images I had used for the first algorithm, in the 4 colors I was using. A convolutional neural network takes in images directly in the form of matrices of pixel values. It then chooses a “filter” matrix and convolves the filter matrix with the image matrices, processes the output, and optimizes the filter to best reproduce the classifications fed in as part of the training data. Again, my main challenge came from the low number of rare supernovae in my sample. I addressed the issue with a different type of oversampling, this time creating additional images by rotating the original images. I still have a long way to go in improving this algorithm, such as by refining the architecture of the neural network, improving my oversampling techniques, and feeding in additional data such as host galaxy properties detected in my previous algorithm, etc. Still, over the summer I was able to achieve 48% accuracy, which is again much better than the 20% random baseline and is already better than my random forest classifier. I believe this algorithm has the potential to achieve significantly higher accuracy with more improvement. <p>
My algorithms will be published in The Astrophysics Journal, a peer-reviewed publication. They will be used in combination with photometric supernova classification algorithms to classify supernovae which do not get spectroscopically classified due to time and resource constraints, to allow for further astrophysics research on these supernovae.
