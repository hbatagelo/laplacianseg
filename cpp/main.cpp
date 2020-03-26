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

#include <QCoreApplication>
#include <QImage>

#include "commandlineparser.h"
#include "seededsegmentationsoft.h"
#include "seededsegmentationhard.h"
#include "seededspsegmentation.h"

int main(int argc, char *argv[])
{
    QCoreApplication app(argc, argv);
    QCoreApplication::setApplicationName("lcseg");
    QCoreApplication::setApplicationVersion("1.0");

    CommandLineParser parser;
    switch (parser.status())
    {
    case CommandLineParser::CommandLineOk:
        break;
    case CommandLineParser::CommandLineError:
        std::cerr << qPrintable(parser.statusMessage())
                  << "\n\n"
                  << qPrintable(parser.helpText());
        return 1;
    case CommandLineParser::CommandLineVersionRequested:
        std::cout << qPrintable(QCoreApplication::applicationName())
                  << qPrintable(QCoreApplication::applicationVersion())
                  << "\n";
        return 0;
    case CommandLineParser::CommandLineHelpRequested:
        parser.showHelp();
        Q_UNREACHABLE();
    }

    // Load input image
    QImage inputImage;
    if (!inputImage.load(parser.inputImagePath()))
    {
        std::cerr << "Cannot load input image file "
                  << qPrintable(parser.inputImagePath()) << ".\n";
        return 0;
    }

    // Load seeds image
    QImage seedsImage;
    if (!seedsImage.load(parser.inputSeedsPath()))
    {
        std::cerr << "Cannot load input seeds image file "
                  << qPrintable(parser.inputSeedsPath()) << ".\n";
        return 0;
    }

    QImage outputImage = QImage(inputImage.width(),
                                inputImage.height(),
                                QImage::Format_ARGB32);

    if (!parser.quietMode())
    {
        std::cout << "Using "
                  << (parser.useHardLabels() ? "hard" : "soft")
                  << " constraints with "
                  << (parser.useSuperpixels() ? "superpixel" : "pixel")
                  << "-based LC.\n";
        std::cout << "Seed colors: "
                  << qPrintable(parser.foregroundSeedColor().name())
                  << " (foreground), "
                  << qPrintable(parser.backgroundSeedColor().name())
                  << " (background).\n";
        if (parser.useSuperpixels())
        {
            std::cout << "Superpixel size " << parser.superpixelSize() << ", ";
            std::cout << "compactness " << parser.compactness() << ".\n";
        }
    }

    // Apply selected segmentation method (LC, LCH, LCSP, LCHSP)
    if (parser.useSuperpixels())
    {
        SeededSPSegmentation seededSeg(&inputImage,
                                       &seedsImage,
                                       parser.superpixelSize(),
                                       parser.compactness());
        seededSeg.Compute(parser.useHardLabels(),
                          parser.foregroundSeedColor(),
                          parser.backgroundSeedColor());
        seededSeg.GetOutput(outputImage, parser.outputAsBinary());
    }
    else
    {
        if (parser.useHardLabels())
        {
            SeededSegmentationHard seededSeg(&inputImage,
                                             &seedsImage,
                                             parser.foregroundSeedColor(),
                                             parser.backgroundSeedColor());
            seededSeg.Compute();
            seededSeg.GetOutput(outputImage, parser.outputAsBinary());
        }
        else
        {
            SeededSegmentationSoft seededSeg(&inputImage,
                                             &seedsImage,
                                             parser.foregroundSeedColor(),
                                             parser.backgroundSeedColor());
            seededSeg.Compute();
            seededSeg.GetOutput(outputImage, parser.outputAsBinary());
        }
    }

    // Save result
    QString name = parser.outputPath();
    if(!outputImage.save(name))
    {
        std::cerr << "Cannot write file " << qPrintable(name) << ".\n";
        return 0;
    }

    if (!parser.quietMode())
    {
        std::cout << "Done! Saved to "
                  << qPrintable(name)
                  << (parser.outputAsBinary() ? " (binary image)" : "")
                  << "\n";
    }
}
