#include "TLorentzVector.h"
#include "TRotation.h"

struct Angles {
  double costh;
  double phi;
};

// NOTE: using TLorentzVector here, although it is possibly way slower than the alternative
Angles calcAnglesHX(const TLorentzVector jpsi, TLorentzVector muPlus)
{
  constexpr double overpi = 1 / M_PI;
  constexpr double pbeam = 6500; // GeV
  constexpr double Mprot = 0.9382720; // GeV
  constexpr double Ebeam = std::sqrt(pbeam*pbeam + Mprot*Mprot);
  TLorentzVector beam1(0, 0, pbeam, Ebeam);
  TLorentzVector beam2(0, 0, -pbeam, Ebeam);

  const auto labToJpsi = -jpsi.BoostVector();
  beam1.Boost(labToJpsi);
  beam2.Boost(labToJpsi);
  muPlus.Boost(labToJpsi);

  const auto beam1Dir = beam1.Vect().Unit(); // in Jpsi frame
  const auto beam2Dir = beam2.Vect().Unit(); // in Jpsi frame
  const auto jpsiDir = jpsi.Vect().Unit(); // in lab frame

  const auto Yaxis = (beam1Dir.Cross(beam2Dir)).Unit() * (jpsi.Rapidity() < 0 ? -1 : 1);

  TRotation rotation;
  rotation.RotateAxes(Yaxis.Cross(jpsiDir), Yaxis, jpsiDir);
  rotation.Invert();

  const auto muPlusRotated = muPlus.Vect().Transform(rotation);

  return Angles{muPlusRotated.CosTheta(), muPlusRotated.Phi() * 180 * overpi};
}
