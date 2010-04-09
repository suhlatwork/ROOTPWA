///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
//
//    This file is part of rootpwa
//
//    rootpwa is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    rootpwa is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rootpwa. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
// File and Version Information:
// $Rev::                             $: revision of last commit
// $Author::                          $: author of last commit
// $Date::                            $: date of last commit
//
// Description:
//      basic test program for amplitude classes
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <fstream>

#include "TVector3.h"
#include "TLorentzRotation.h"

#include "Vec.h"
#include "lorentz.h"
#include "keyfile.h"
// #include "particle.h"
#include "event.h"
#include "wave.h"

#include "utilities.h"
#include "particleDataTable.h"
#include "diffractiveDissVertex.h"
#include "diffractiveDissVertex2.h"
#include "massDependence.h"
#include "isobarHelicityAmplitude.h"
#include "isobarHelicityAmplitude2.h"


extern particleDataTable PDGtable;
extern wave              gWave;


using namespace std;
using namespace rpwa;


int
main(int argc, char** argv)
{

  rpwa::particleDataTable& pdt = rpwa::particleDataTable::instance();
  pdt.readFile();

  rpwa::isobarDecayVertex::setDebug(true);
  rpwa::isobarDecayTopology::setDebug(true);
  rpwa::isobarHelicityAmplitude::setDebug(true);

  rpwa::isobarDecayVertex2::setDebug(true);
  rpwa::isobarDecayTopology2::setDebug(true);
  rpwa::isobarHelicityAmplitude2::setDebug(true);

  if (0) {
    {
      fourVec  p(2, threeVec(0.5, 0.75, 1));
      threeVec n = threeVec(0, 0, 1) / p.V();
      cout << "before = " << n << "    " << p << endl;
      rotation         R;
      lorentzTransform L1;
      L1.set(R.set(n.phi(), n.theta() - M_PI / 2, -M_PI / 2));
      n *= R;
      p *= L1;
      cout << "L1 -> " << n << "    " << p << endl;
      lorentzTransform L2;
      L2.set(R.set(0, signof(p.x()) * p.theta(), 0));
      p *= L2;
      cout << "L2 -> " << p << endl;
      lorentzTransform L3;
      L3.set(p);
      p *= L3;
      cout << "L3 -> " << p << endl;

      matrix<double> X(4, 4);
      X = L3 * (L2 * L1);
      lorentzTransform L(X);
      p = fourVec(2, threeVec(0.5, 0.75, 1));
      p *= L;
      cout << "L -> " << p << endl;
    }

    {
      TLorentzVector p(0.5, 0.75, 1, 2);
      TVector3       n = TVector3(0, 0, 1).Cross(p.Vect());
      TRotation R1;
      R1.RotateZ(-n.Phi());
      R1.RotateY(piHalf - n.Theta());
      R1.RotateZ(piHalf);
      n *= R1;
      p *= R1;
      cout << "R1 -> " << n << "    " << p << endl;
      // rotate about yHfAxis so that daughter momentum is along z-axis
      TRotation R2;
      R2.RotateY(-signum(p.X()) * p.Theta());
      p *= R2;
      cout << "R2 -> " << p << endl;
      // boost into daughter RF
      TLorentzRotation L3(-p.BoostVector());
      cout << "L3 -> " << L3 * p << endl;

      R1.Transform(R2);
      TLorentzRotation L(R1);
      L.Boost(-p.BoostVector());
      p = TLorentzVector(0.5, 0.75, 1, 2);
      p *= L;
      cout << "L -> " << p << endl;
    }

    {
      rpwa::particle X("X");
      TLorentzVector p(0.5, 0.75, 1, 2);
      X.setLzVec(p);
      rpwa::isobarHelicityAmplitude amp;
      TLorentzRotation L = amp.hfTransform(X.lzVec());
      cout << "!!! L -> " << L * p << endl;
    }
  }

  if (0) {
    TLorentzVector beam(1,   0.5,  180, 182);
    TLorentzVector X   (0.5, 0.75, 1,   3);
    rpwa::isobarHelicityAmplitude amp;
    TLorentzRotation L = amp.gjTransform(beam, X);
    cout << "!!! L -> " << L * X << endl;
  }

  if (0) {
    // define final state particles
    rpwa::particle pi0("pi-");
    rpwa::particle pi1("pi+");
    rpwa::particle pi2("pi-");
    rpwa::particle pi3("pi+");
    rpwa::particle pi4("pi-");
    // define isobars
    rpwa::particle sigma("sigma");
    rpwa::particle a1   ("a1(1269)+");
    rpwa::particle f1   ("f1(1285)");
    // define X-system
    //                     I   G  2J  P   C  2M
    rpwa::particle X("X-", 2, -1, 4, +1, +1, 2);
    f1.setSpinProj(-2);
    a1.setSpinProj(-2);
    // define production vertex
    rpwa::particle              beam("pi-");
    rpwa::diffractiveDissVertex prodVert(beam, X);
    // define vertices
    rpwa::isobarDecayVertex vert0(X,     pi4, f1,    2, 2);
    rpwa::isobarDecayVertex vert1(f1,    pi2, a1,    2, 2);
    rpwa::isobarDecayVertex vert2(a1,    pi3, sigma, 2, 0);
    rpwa::isobarDecayVertex vert3(sigma, pi0, pi1,   0, 0);
    // set Lorentz vectors
    beam.setLzVec(TLorentzVector(0.104385398, 0.0132061851, 189.987978, 189.988058));
    pi0.setLzVec(TLorentzVector(-0.0761465106, -0.116917817, 5.89514709, 5.89844947));
    pi1.setLzVec(TLorentzVector(-0.0244305532, -0.106013023, 30.6551865, 30.6556973));
    pi2.setLzVec(TLorentzVector(0.000287952441, 0.10263611, 3.95724077, 3.96103114));
    pi3.setLzVec(TLorentzVector(0.0299586212, 0.176440177, 115.703054, 115.703277));
    pi4.setLzVec(TLorentzVector(0.176323963, -0.0985753246, 30.9972271, 30.9981995));
    // build graph
    vector<rpwa::particle*> fsParticles;
    fsParticles.push_back(&pi0);
    fsParticles.push_back(&pi1);
    fsParticles.push_back(&pi2);
    fsParticles.push_back(&pi3);
    fsParticles.push_back(&pi4);
    vector<rpwa::isobarDecayVertex*> decayVertices;
    decayVertices.push_back(&vert3);
    decayVertices.push_back(&vert1);
    decayVertices.push_back(&vert2);
    decayVertices.push_back(&vert0);
    rpwa::isobarDecayTopology     topo(fsParticles, prodVert, decayVertices);
    rpwa::isobarHelicityAmplitude amp(topo);

    // cout << "!!! before: "<< endl;
    // for (unsigned int i = 0; i < topo.nmbVertices(); ++i)
    //   cout << *topo.isobarDecayVertices()[i];
    // amp.transformDaughters();
    // cout << "!!! after: " << endl;
    // for (unsigned int i = 0; i < topo.nmbVertices(); ++i)
    //   cout << *topo.isobarDecayVertices()[i];
    // for (unsigned int i = 0; i < topo.nmbVertices(); ++i) {
    //   complex<double> decayAmp = amp.twoBodyDecayAmplitude(*topo.isobarDecayVertices()[i], (i == 0));
    //   cout << topo.isobarDecayVertices()[i]->mother().name() << " decay amplitude = "
    // 	   << decayAmp << endl;
    // }

    // complex<double> decayAmp = amp.twoBodyDecayAmplitudeSum(*topo.isobarDecayVertices()[0], true);
    complex<double> decayAmp = amp.amplitude();
    cout << "!!!< decay amplitude = " << decayAmp << endl;
    // cout << beam << endl;
    // for (unsigned int i = 0; i < fsParticles.size(); ++i)
    //   cout << *fsParticles[i] << endl;

  }

  if (1) {
    // define final state particles
    particlePtr pi0 = createParticle("pi-");
    particlePtr pi1 = createParticle("pi+");
    particlePtr pi2 = createParticle("pi-");
    particlePtr pi3 = createParticle("pi+");
    particlePtr pi4 = createParticle("pi-");
    // define isobars
    particlePtr sigma = createParticle("sigma");
    particlePtr a1    = createParticle("a1(1269)+");
    particlePtr f1    = createParticle("f1(1285)");
    // define X-system
    //                     I   G  2J  P   C  2M
    particlePtr X = createParticle("X-", 2, -1, 4, +1, +1, 2);
    f1->setSpinProj(-2);
    a1->setSpinProj(-2);
    // define production vertex
    particlePtr              beam     = createParticle("pi-");
    diffractiveDissVertexPtr prodVert = createDiffractiveDissVertex(beam, X);
    // define vertices
    massDependencePtr    massDep = createRelativisticBreitWigner();
    isobarDecayVertexPtr vert0   = createIsobarDecayVertex(X,     pi4, f1,    2, 2);
    isobarDecayVertexPtr vert1   = createIsobarDecayVertex(f1,    pi2, a1,    2, 2, massDep);
    isobarDecayVertexPtr vert2   = createIsobarDecayVertex(a1,    pi3, sigma, 2, 0, massDep);
    isobarDecayVertexPtr vert3   = createIsobarDecayVertex(sigma, pi0, pi1,   0, 0, massDep);
    // set Lorentz vectors
    beam->setLzVec(TLorentzVector(0.104385398, 0.0132061851, 189.987978, 189.988058));
    pi0->setLzVec(TLorentzVector(-0.0761465106, -0.116917817, 5.89514709, 5.89844947));
    pi1->setLzVec(TLorentzVector(-0.0244305532, -0.106013023, 30.6551865, 30.6556973));
    pi2->setLzVec(TLorentzVector(0.000287952441, 0.10263611, 3.95724077, 3.96103114));
    pi3->setLzVec(TLorentzVector(0.0299586212, 0.176440177, 115.703054, 115.703277));
    pi4->setLzVec(TLorentzVector(0.176323963, -0.0985753246, 30.9972271, 30.9981995));
    // build graph
    vector<isobarDecayVertexPtr> decayVertices;
    decayVertices.push_back(vert3);
    decayVertices.push_back(vert1);
    decayVertices.push_back(vert2);
    decayVertices.push_back(vert0);
    vector<particlePtr> fsParticles;
    fsParticles.push_back(pi0);
    fsParticles.push_back(pi1);
    fsParticles.push_back(pi2);
    fsParticles.push_back(pi3);
    fsParticles.push_back(pi4);
    isobarDecayTopology2     topo(prodVert, decayVertices, fsParticles);
    isobarHelicityAmplitude2 amp(topo);
    complex<double>          decayAmp = amp.amplitude();
    cout << "!!!< decay amplitude = " << decayAmp << endl;

    if (1) {  // compare to PWA2000
      PDGtable.initialize();
      event    ev;
      ifstream eventData("testEvents.evt");
      ev.setIOVersion(1);
      if (!(eventData >> ev).eof()) {
	keyfile key;
	key.open("1-2-+1+pi-_11_f11285=pi-_11_a11269=pi+_1_sigma.key");
	complex<double> pwa2000amp;
	key.run(ev, pwa2000amp, true);
	key.rewind();
	key.close();
	cout << "!!! PWA2000 amplitude = " << pwa2000amp << endl;
	cout << "!!! my amplitude = " << decayAmp << " vs. PWA2000 amplitude = " << pwa2000amp << ", "
	     << "delta = " << decayAmp - pwa2000amp << endl;
      }
    }
  }

}
