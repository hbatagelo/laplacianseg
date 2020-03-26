/****************************************************************************
**
** Copyright (C) 2020 Harlen Batagelo <hbatagelo@gmail.com> and
**                    João Paulo Gois <jpgois@gmail.com>.
**
** This file is part of the implementation of the paper
** 'Laplacian Coordinates: Theory and Methods for Seeded Image Segmentation'
** by Wallace Casaca, João Paulo Gois, Harlen Batagelo, Gabriel Taubin and
** Luis Gustavo Nonato.
**
** This program is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with this program. If not, see <http://www.gnu.org/licenses/>.
**
****************************************************************************/

#include "opencv2/opencv.hpp"
#include "seededsegmentationhard.h"

SeededSegmentationHard::SeededSegmentationHard(QImage *inputImage,
                                               QImage *seedsImage,
                                               QColor fgColor,
                                               QColor bgColor,
                                               double beta)
    : m_image(inputImage),
      m_beta(beta)
{
    m_cols = inputImage->width();
    m_rows = inputImage->height();
    m_pixelCount = m_rows * m_cols;

    m_inputImage[0] = Eigen::MatrixXd(m_rows, m_cols);
    m_inputImage[1] = Eigen::MatrixXd(m_rows, m_cols);
    m_inputImage[2] = Eigen::MatrixXd(m_rows, m_cols);
    m_seedsMatrix = Eigen::MatrixXi(m_rows, m_cols);

    // Split the different layers of the image and store them in
    // three Eigen matrices
    auto stInput = reinterpret_cast<QRgb *>(inputImage->bits());
    auto stSeeds = reinterpret_cast<QRgb *>(seedsImage->bits());
    for (auto i = 0; i < m_rows; ++i)
    {
        for (auto j = 0; j < m_cols; ++j)
        {
            auto offset = i * m_cols + j;

            auto inputColor = QColor(stInput[offset]);
            m_inputImage[0](i, j) = inputColor.redF();
            m_inputImage[1](i, j) = inputColor.greenF();
            m_inputImage[2](i, j) = inputColor.blueF();

            auto seedColor = QColor(stSeeds[offset]);
            if (seedColor == fgColor)
            {
                m_seedsMatrix(i, j) = 1;
                m_idxBackground.push_back(offset);
            }
            else if (seedColor == bgColor)
            {
                m_seedsMatrix(i, j) = 0;
                m_idxForeground.push_back(offset);
            }
            else
            {
                m_idxUnknown.push_back(offset);
                m_seedsMatrix(i, j) = -1;
            }
        }
    }
}

// Initialize the sparse matrices with the number of expected values in
// each row and instantiate the variables which will help us to compute
// the graph weight and the sum of the rows
void SeededSegmentationHard::SetUpGraphWeightAndSum()
{
    const double epsilon = 10e-6;

    m_sparseWeights = Eigen::SparseMatrix<double>(m_pixelCount, m_pixelCount);
    m_sparseWeights.reserve(Eigen::VectorXi::Constant(m_pixelCount, 9));

    m_sparseD = Eigen::SparseMatrix<double>(m_pixelCount, m_pixelCount);
    m_sparseD.reserve(Eigen::VectorXi::Constant(m_pixelCount, 1));

    // Compute weights of every pixels with their neighbors and uptate the sum
    // of the row weights at every iterations
    for (int i = 0; i < m_pixelCount; i++)
    {
        // Get row and column of the pixel
        double summ = 0;
        int rowPi = i / m_cols;
        int colPi = i % m_cols;

        for (int j = 0; j < 9; j++)
        {
            // Get row an column of the neighbor
            int rowPj = rowPi + j / 3 - 1;
            int colPj = colPi + j % 3 - 1;

            int indexJ = (rowPj * m_cols) + colPj;

            // If neighbor is out of the image bounds or is corresponding to
            // the same pixel, do nothing
            if (!(rowPj < 0 || rowPj >= m_rows ||
                  colPj < 0 || colPj >= m_cols ||
                  (rowPj == rowPi && colPj == colPi)))
            {
                // Compute difference of the values of
                // each layer between the two pixels
                Eigen::VectorXd VecDiff(3);
                for (auto k = 0; k < 3; ++k)
                {
                    VecDiff(k) = m_inputImage[k](rowPi, colPi) -
                            m_inputImage[k](rowPj, colPj);
                }
                // Calculate infinity norm of the differences
                double dist = VecDiff.lpNorm<Eigen::Infinity>();
                // Compute the weight
                double wij = std::exp(-m_beta * dist * dist);
                wij = std::round(wij * 1e6) * 1e-6;
                wij += epsilon;
                m_sparseWeights.insert(i, indexJ) = wij;

                summ += wij;
            }
        }
        // Insert the sum in the diagonal matrix
        m_sparseD.insert(i, i) = summ;
    }
}

// Compute difference between the row sum of the weights and the weights
void SeededSegmentationHard::DiffBetweenSumAndWeights()
{
    m_sparseL = Eigen::SparseMatrix<double>(m_pixelCount, m_pixelCount);
    m_sparseL.reserve(Eigen::VectorXi::Constant(m_pixelCount, 9));
    m_sparseL = m_sparseD - m_sparseWeights;
}

bool SeededSegmentationHard::Compute()
{
    SetUpGraphWeightAndSum();
    DiffBetweenSumAndWeights();
    SolveHardLinearSystem();

    return true;
}

void SeededSegmentationHard::SolveHardLinearSystem()
{
    std::vector<int> seededidx;
    seededidx.insert(seededidx.end(), m_idxBackground.begin(), m_idxBackground.end());
    seededidx.insert(seededidx.end(), m_idxForeground.begin(), m_idxForeground.end());

    std::vector<int> vectorpermutation;
    vectorpermutation.insert(vectorpermutation.end(), seededidx.begin(), seededidx.end());
    vectorpermutation.insert(vectorpermutation.end(), m_idxUnknown.begin(), m_idxUnknown.end());

    Eigen::VectorXi permutation;
    permutation = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(
                vectorpermutation.data(),
                long(vectorpermutation.size()));
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic,int> perm(permutation);

    Eigen::SparseMatrix<double> spLtwisted;
    spLtwisted = m_sparseL.twistedBy(perm.transpose());

    auto b_rows = long(seededidx.size());

    Eigen::SparseMatrix<double> sparseLs;
    sparseLs = spLtwisted.topLeftCorner(b_rows, b_rows);

    auto b_cols = m_sparseL.cols() - b_rows;
    Eigen::SparseMatrix<double> sparseB;
    sparseB= spLtwisted.block(0, b_rows, b_rows, b_cols);
    Eigen::SparseMatrix<double> sparseBt = sparseB.transpose();

    Eigen::SparseMatrix<double> sparseLu;
    sparseLu = spLtwisted.block(b_rows, b_rows, b_cols, b_cols);

    Eigen::VectorXd xs(seededidx.size());
    xs.setConstant(1);
    xs.segment(0, long(m_idxBackground.size())).setConstant(0);

    Eigen::SparseMatrix<double> sparseLHard;
    sparseLHard = sparseBt * sparseB + sparseLu * sparseLu;
    Eigen::VectorXd bHard = -(sparseLu * sparseBt + sparseBt * sparseLs) * xs;

    Eigen::CholmodDecomposition<Eigen::SparseMatrix<double>, Eigen::Lower> sparseSolver;

    sparseSolver.setMode(Eigen::CholmodSupernodalLLt);
    sparseSolver.compute(sparseLHard);
    m_solutionVector = sparseSolver.solve(bHard);
}

void SeededSegmentationHard::GetOutput(QImage &outputImage, bool asBinary)
{
    int width = m_image->width();
    int height = m_image->height();

    bool otsu = true;
    // If otsuVariant is true, the threshold is
    // given by Isol > (1-otsu) instead of Isol > otsu
    bool otsuVariant = true;

    cv::Mat imBw(height, width, CV_8UC1);
    cv::Mat imGray(height, width, CV_8UC1);

    // Mount solution vector as an OpenCV matrix
    for (auto i = 0; i < m_rows; ++i)
    {
        for (auto j = 0; j < m_cols; ++j)
        {
            double value = 0;
            int seedType = m_seedsMatrix(i, j);
            if(seedType == 1)
                value = 255;

            auto cpt = cv::Point(j, i);
            imGray.at<uchar>(cpt) = uchar(value);
        }
    }

    for(int i = 0; i < m_solutionVector.size(); i++)
    {
        double value = 1 - m_solutionVector.coeff(i);

        value *= 255.0;
        value = std::clamp(value, 0.0, 255.0);
        if (std::isnan(value))
            value = 0;

        int offset = m_idxUnknown[ulong(i)];
        int rowPi = offset / m_cols;
        int colPi = offset % m_cols;

        auto cpt = cv::Point(colPi, rowPi);
        imGray.at<uchar>(cpt) = uchar(value);
    }

    // Threshold with Otsu
    if (otsu)
    {
        if (otsuVariant)
        {
            double otsuThreshold = cv::threshold(
                imGray, imBw, 0, 255, CV_THRESH_BINARY | CV_THRESH_OTSU);
            if (almostEquals(otsuThreshold, 0.0))
                otsuThreshold = 255;

            cv::threshold(imGray, imBw, 255 - otsuThreshold, 255, CV_THRESH_BINARY);
        }
        else
        {
            cv::threshold(imGray, imBw, 0, 255, CV_THRESH_BINARY | CV_THRESH_OTSU);
        }
    }
    else
    {
        cv::threshold(imGray, imBw, 127, 255, CV_THRESH_BINARY);
    }

    for (auto i = 0; i < height; ++i)
    {
        for (auto j = 0; j < width; ++j)
        {
            auto cpt = cv::Point(j, i);
            auto pt = QPoint(j, i);

            if(imBw.at<uchar>(cpt) == 255)
            {
                if (asBinary)
                    outputImage.setPixelColor(pt, QColor(255,255,255,255));
                else
                    outputImage.setPixelColor(pt, m_image->pixel(j,i));
            }
            else
            {
                if (asBinary)
                    outputImage.setPixelColor(pt, QColor(0,0,0,255));
                else
                    outputImage.setPixelColor(pt, QColor(0,0,0,0));
            }
        }
    }
}
