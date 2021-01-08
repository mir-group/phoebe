#include "parser.h"

std::tuple<Crystal, PhononH0> Parser::parsePhHarmonic(Context &context) {
  std::string fileName = context.getPhD2FileName();
  // TODO -- this feels like a problem waiting to happen. For example, 
  // what if someone named some part of their filepath ".fc"?

  // check if this is a qe fc file
  if (fileName.find(".fc") != std::string::npos) {
    return QEParser::parsePhHarmonic(context);
  }
  else {
    return PhonopyParser::parsePhHarmonic(context);
  }
}

std::tuple<Crystal, ElectronH0Fourier> Parser::parseElHarmonicFourier(
        Context &context) {
  return QEParser::parseElHarmonicFourier(context);
}

std::tuple<Crystal, ElectronH0Wannier> Parser::parseElHarmonicWannier(
            Context &context, Crystal *inCrystal) {
  return QEParser::parseElHarmonicWannier(context, inCrystal);
}
