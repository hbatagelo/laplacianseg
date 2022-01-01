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

#include "commandlineparser.h"

CommandLineParser::CommandLineParser() {
  m_parser.setApplicationDescription(
      "Seeded image segmentation using Laplacian Coordinates.\n"
      "(c) 2020-2022 Harlen Batagelo, João Paulo Gois\n\n"
      "From the paper \"Casaca et al., Laplacian Coordinates: Theory and "
      "Methods for Seeded Image Segmentation, DOI "
      "10.1109/TPAMI.2020.2974475.\"");

  m_parser.setSingleDashWordOptionMode(QCommandLineParser::ParseAsLongOptions);

  m_parser.addPositionalArgument("input", "Input image file.");
  m_parser.addPositionalArgument("seeds", "Input seeds image file.");
  m_parser.addPositionalArgument("output", "Output image file.");

  QCommandLineOption const foregroundColorOption{
      "fg", "Sets color of foreground seeds (default is '#ff0000').", "color",
      "#ff0000"};
  m_parser.addOption(foregroundColorOption);

  QCommandLineOption const backgroundColorOption{
      "bg", "Sets color of background seeds (default is '#0000ff').", "color",
      "#0000ff"};
  m_parser.addOption(backgroundColorOption);

  QCommandLineOption const hardOption{
      "hard", "Uses seeds as hard labeling constraints."};
  m_parser.addOption(hardOption);

  QCommandLineOption const superpixelOption{"superpixel",
                                            "Uses SLIC superpixels."};
  m_parser.addOption(superpixelOption);

  QCommandLineOption const superpixelSizeOption{
      "size", "Sets superpixel size (default is 100).", "value", "100"};
  m_parser.addOption(superpixelSizeOption);

  QCommandLineOption const superpixelCompactnessOption{
      "compactness", "Sets superpixel compactness (default is 10.0).", "value",
      "10.0"};
  m_parser.addOption(superpixelCompactnessOption);

  QCommandLineOption const binaryOption{{"b", "binary"},
                                        "Writes output as a binary image."};
  m_parser.addOption(binaryOption);

  QCommandLineOption const quietOption{{"q", "quiet"}, "Runs in silent mode."};
  m_parser.addOption(quietOption);

  QCommandLineOption const versionOption{m_parser.addVersionOption()};

  if (!m_parser.parse(QCoreApplication::arguments())) {
    m_statusMessage = m_parser.errorText();
    m_status = CommandLineError;
    return;
  }

  if (m_parser.isSet(versionOption)) {
    m_status = CommandLineVersionRequested;
    return;
  }

  QStringList const positionalArguments{m_parser.positionalArguments()};
  if (positionalArguments.isEmpty()) {
    m_statusMessage = "Argument 'input' missing.";
    m_status = CommandLineError;
    return;
  }

  if (positionalArguments.size() == 1) {
    m_statusMessage = "Argument 'seeds' missing.";
    m_status = CommandLineError;
    return;
  }

  if (positionalArguments.size() == 2) {
    m_statusMessage = "Argument 'output' missing.";
    m_status = CommandLineError;
    return;
  }

  QStringList const arguments{m_parser.positionalArguments()};
  m_inputImagePath = arguments.at(0);
  m_inputSeedsPath = arguments.at(1);
  m_outputPath = arguments.at(2);
  m_useHardLabels = m_parser.isSet(hardOption);
  m_useSuperpixels = m_parser.isSet(superpixelOption);
  m_outputAsBinary = m_parser.isSet(binaryOption);
  m_quietMode = m_parser.isSet(quietOption);
  m_superpixelSize = m_parser.value(superpixelSizeOption).toInt();
  if (m_superpixelSize < 1) {
    m_statusMessage = "Invalid superpixel size.";
    m_status = CommandLineError;
    return;
  }
  m_compactness = m_parser.value(superpixelCompactnessOption).toDouble();
  if (m_compactness < 1.0) {
    m_statusMessage = "Invalid compactness.";
    m_status = CommandLineError;
    return;
  }

  QString const foregroundColorName{m_parser.value(foregroundColorOption)};
  if (QColor::isValidColor(foregroundColorName)) {
    m_foregroundSeedColor = QColor(foregroundColorName);
  } else {
    m_statusMessage = "Invalid foreground seed color.";
    m_status = CommandLineError;
    return;
  }

  QString const backgroundColorName{m_parser.value(backgroundColorOption)};
  if (QColor::isValidColor(backgroundColorName)) {
    m_backgroundSeedColor = QColor(backgroundColorName);
  } else {
    m_statusMessage = "Invalid background seed color.";
    m_status = CommandLineError;
    return;
  }
}
