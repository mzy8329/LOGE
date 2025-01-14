#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <opencv2/opencv.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/features2d.hpp>

#include "LOGE.h"

std::vector<cv::Mat> readImagesFromFolder(const std::string& folderPath) {
    std::vector<cv::Mat> images;
    for (int i = 0; ; ++i)
    {
        std::string imagePath = folderPath + "/" + std::to_string(i) + ".jpg";
        cv::Mat img = cv::imread(imagePath);
        if (img.empty())
        {
            break;
        }
        images.push_back(img);
    }
    return images;
}



int main()
{
    std::string folderPath = "../data/imgs";
    std::vector<cv::Mat> images = readImagesFromFolder(folderPath);

    cv::Ptr<cv::ORB> detector = cv::ORB::create(300);
    for (const auto& image : images)
    {
        cv::Mat gray;
        cv::cvtColor(image, gray, cv::COLOR_BGR2GRAY);
        std::vector<cv::KeyPoint> keyPoints;
        detector->detect(gray, keyPoints);


        cv::Mat img_with_kp;
        cv::drawKeypoints(image, keyPoints, img_with_kp, cv::Scalar::all(-1), cv::DrawMatchesFlags::DEFAULT);

        cv::imshow("Image", img_with_kp);
        cv::waitKey(0);
    }

    return 0;
}