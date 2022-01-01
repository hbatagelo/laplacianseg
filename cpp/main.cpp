/****************************************************************************
**
** Copyright (C) 2020-2022 Harlen Batagelo <hbatagelo@gmail.com>,
**                         João Paulo Gois <jpgois@gmail.com>.
**
** This file is part of the implementation of the paper
** 'Laplacian Coordinates: Theory and Methods for Seeded Image Segmentation'
** by Wallace Casaca, João Paulo Gois, Harlen Batagelo, Gabriel Taubin and
** Luis Gustavo Nonato. DOI 10.1109/TPAMI.2020.2974475.
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
#include <opencv2/core/utils/logger.hpp>

#include "commandlineparser.h"
#include "seededsegmentationhard.h"
#include "seededsegmentationsoft.h"
#include "seededspsegmentation.h"

int main(int argc, char *argv[]) {
  QCoreApplication app{argc, argv};
  QCoreApplication::setApplicationName("lcseg");
  QCoreApplication::setApplicationVersion("1.0.1");

  cv::utils::logging::setLogLevel(
      cv::utils::logging::LogLevel::LOG_LEVEL_SILENT);

  CommandLineParser parser;
  switch (parser.status()) {
  case CommandLineParser::CommandLineOk:
    break;
  case CommandLineParser::CommandLineError:
    std::cerr << qPrintable(parser.statusMessage()) << "\n\n"
              << qPrintable(parser.helpText());
    return 1;
  case CommandLineParser::CommandLineVersionRequested:
    std::cout << qPrintable(QCoreApplication::applicationName())
              << qPrintable(QCoreApplication::applicationVersion()) << "\n";
    return 1;
  case CommandLineParser::CommandLineHelpRequested:
    parser.showHelp();
    Q_UNREACHABLE();
  }

  // Load input image
  QImage input;
  QString inputImagePath{parser.inputImagePath()};
  if (!input.load(inputImagePath)) {
    std::cerr << "Cannot load input image file "
              << qPrintable(parser.inputImagePath()) << ".\n";
    return 1;
  }

  // Load seeds image
  QImage seeds;
  if (!seeds.load(parser.inputSeedsPath())) {
    std::cerr << "Cannot load input seeds image file "
              << qPrintable(parser.inputSeedsPath()) << ".\n";
    return 1;
  }

  if (input.width() != seeds.width() || input.height() != seeds.height()) {
    std::cerr << "Input image and seeds image must have the same size.\n";
    return 1;
  }

  if (!parser.quietMode()) {
    std::cout << "Using " << (parser.useHardLabels() ? "hard" : "soft")
              << " constraints with "
              << (parser.useSuperpixels() ? "superpixel" : "pixel")
              << "-based LC.\n";
    std::cout << "Seed colors: "
              << qPrintable(parser.foregroundSeedColor().name())
              << " (foreground), "
              << qPrintable(parser.backgroundSeedColor().name())
              << " (background).\n";
    if (parser.useSuperpixels()) {
      std::cout << "Superpixel size " << parser.superpixelSize() << ", ";
      std::cout << "compactness " << parser.compactness() << ".\n";
    }
  }

  QImage output{input.width(), input.height(), QImage::Format_ARGB32};

  try {
    // Apply selected segmentation method (LC, LCH, LCSP, LCHSP)
    if (parser.useSuperpixels()) {
      lc::SeededSPSegmentation seededSeg(
          &input, &seeds, parser.superpixelSize(), parser.compactness());
      seededSeg.Compute(parser.useHardLabels(), parser.foregroundSeedColor(),
                        parser.backgroundSeedColor());
      seededSeg.GetOutput(output, parser.outputAsBinary());
    } else {
      if (parser.useHardLabels()) {
        lc::SeededSegmentationHard seededSeg(&input, &seeds,
                                             parser.foregroundSeedColor(),
                                             parser.backgroundSeedColor());
        seededSeg.Compute();
        seededSeg.GetOutput(output, parser.outputAsBinary());
      } else {
        lc::SeededSegmentationSoft seededSeg(&input, &seeds,
                                             parser.foregroundSeedColor(),
                                             parser.backgroundSeedColor());
        seededSeg.Compute();
        seededSeg.GetOutput(output, parser.outputAsBinary());
      }
    }
  } catch (std::bad_alloc const &e) {
    std::cerr << e.what() << '\n';
  } catch (...) {
    std::cerr << "An exception has occurred.\n";
    return 1;
  }

  // Save result
  QString const name{parser.outputPath()};
  if (!output.save(name)) {
    std::cerr << "Cannot write file " << qPrintable(name) << ".\n";
    return 1;
  }

  if (!parser.quietMode()) {
    std::cout << "Done! Saved to " << qPrintable(name)
              << (parser.outputAsBinary() ? " (as binary image)" : "") << ".\n";
  }
}
