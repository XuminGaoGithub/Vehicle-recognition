#A fast and effective vehicle recognition algorithm based on HOG and Random Forest(July 2019):

Authur：Xumin Gao, Jing Zhao(Tom). (Qingdao Smart Ground Vehicle Intelligent Technology Co. Ltd)

E-mail: comin15071460998@gmail.com

# Requirements

- Cmake 3.1

- boost

- libxml:

sudo apt-get install libxml2-dev libxslt1-dev

sudo su

ln -s /usr/include/libxml2/libxml   /usr/include/libxml

- libconfig libconfig++:

sudo apt-get install libconfig-dev libconfig++-dev libconfig-dev

sudo apt-get install libconfig++-dev libconfig-dev

sudo apt install libgl-dev libglu-dev libglib2.0-dev libsm-dev libxrender-dev libfontconfig1-dev libxext-dev


# Usage

###./bin or RandomForest1:

Executable file generated by cmake.

###./pic:

We convert the original image into a grayscale image and compress the pixel value to 0-1 for storage. This can speed up the loading of data.

###./label:

The label of the corrosponding type of vehicles for ./pic 

###./data:
We converted the features extracted from the vehicle picture into HOG features(This part requires an artificial label), and the vehicle in each picture was extracted form four sub-images.

For example:

Featureout.data:the HOG features from bumpers

Featureoutlabel:label corresponding the type of vehicles


###./include

The library of Random Forest algorithm.

###RandomForest.cpp

Main file(.cpp) in this project

###./src

Other .cpp file in this project

###./xml

Models of the different part of vehicles after training

###slidewindow.conf

Although we have done the label, in order to expand the amount of training data, we use the fine-tuning sliding window to generate more training data on the basis of the original label.


###Run

cd build
cmake .
make

cd ..
./RandomForest1




#Abstract

We developed an algorithm to classify different categories of vehicles. The vehicle was classified by a random forest model, pre-trained by HOG features extracted from different parts in the vehicle images. The accuracy reached 91.26%. Moreover, the classification only needed 0.012s on the CPU for one image with 1024x1024 size. The prediction time was greatly shortened compared with the convolutional neural network. 

After the vehicle was classified, we reconstructed the 3D model of the vehicle by extracting feature points and the depth information of the corresponding category of vehicles.

***We have opened source the code of vehicle recognition. But we still cannot open source the code of vehicle 3D modeling due to confidentiality agreement.

![Image text](https://github.com/XuminGaoGithub/Vehicle-recognition/blob/main/1.png) 




