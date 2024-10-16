#include "uego.h"
#include "time.h"
//#include "internal.h"
#include <stdlib.h>
#include <unistd.h>
#include "Tanimoto.h"
#include <iomanip>

#include "openeye.h"

#include "oezap.h"
#include "oechem.h"
#include "oesystem.h"
#include "oeplatform.h"

#include <string>
#include <fstream>
#include <streambuf>

using namespace OEPB;
using namespace OEChem;
using namespace OESystem;
using namespace OEPlatform;

////////////////////////////////////////////////////////////
// $Id: testval.cc,v 2.6 1998/03/29 10:40:32 jelasity Exp $
// testval.cc
// definition of NDimRealElement::Value()
// objective function library for n dimensional real spaces
// this version was originaly created for experimental tests
////////////////////////////////////////////////////////////
// modification history:
////////////////////////////////////////////////////////////

//-------------------------------------------------------------------------
#include "VolumeOverlap.h"
#include "MoveAndRotate.h"
#include "Tanimoto.h"
#include <iostream>

string GetStdoutFromCommand2(string cmd) {

	string data;
	FILE * stream;
	const int max_buffer = 256;
	char buffer[max_buffer];
	cmd.append(" 2>&1");

	stream = popen(cmd.c_str(), "r");
	if (stream) {
		while (!feof(stream))
			if (fgets(buffer, max_buffer, stream) != NULL)
				data.append(buffer);
		pclose(stream);
	}
	return data;
}

double NDimMolElement::Value() {

	if (x[1] == x[4] && x[2] == x[5] && x[3] == x[6]) {
		return 0;
	}

	double * newAtomsRotated =
			&(INI.getMolVariable()->getAtomsXYZ()[INI.getMolVariable()->atoms.size()
					* 3]);
	MoveAndRotate::RotateMolAccording1Axis(INI.getMolVariable()->getAtomsXYZ(),
	INI.getMolVariable()->atoms.size(), x[0], x[1], x[2], x[3], x[4], x[5],
			x[6], newAtomsRotated);

	MoveAndRotate::MolToNewPosition(newAtomsRotated,
	INI.getMolVariable()->atoms.size(), x[7], x[8], x[9]);

	////////////RESTRICCION SHAPE SIMILARITY
	double VAB = VolumeOverlap::overlapWEGA(
	INI.getMolQuery()->getAtomsXYZ(),
	INI.getMolQuery()->atoms.size(), INI.getMolQuery()->getWeightAtoms(),
	INI.getMolQuery()->getRadiusAtoms(), newAtomsRotated,
	INI.getMolVariable()->atoms.size(),
	INI.getMolVariable()->getWeightAtoms(),
	INI.getMolVariable()->getRadiusAtoms(),
	INI.getSameVanDerWaalsRadius());
	//cout << VAB << endl;
	//cout << INI.getMolQuery()->tanimoto << endl;
	//cout << INI.getMolVariable()->tanimoto << endl;
	double result2 = Tanimoto::calculateTanimotoGeneric(
	INI.getMolQuery()->tanimoto,
	INI.getMolVariable()->tanimoto, VAB);
	//cout << result2 << "\n";

	if (result2 < 0.0001) {

		//cout << "------------------------------------------\n";
		return 0;
	}

	///////////////////
	///////////////////////////////////
	INI.getMoleculeTemp()->updateAtomCoordinates(newAtomsRotated);

	OEGraphMol refmol, fitmol;

	oemolistream ifs;
	ifs.SetFormat(OEFormat::MOL2);
	ifs.openstring(INI.getMolQueryToString());

	OEReadMolecule(ifs, refmol);

	OEMMFFAtomTypes(refmol);
	OEMMFF94PartialCharges(refmol);
	OEAssignBondiVdWRadii(refmol);

	OEET et;
	et.SetRefMol(refmol);

	float tanimoto;

	ifs.openstring(WriteMolecule::MolToString(INI.getMoleculeTemp()));
	OEReadMolecule(ifs, fitmol);
		OEMMFFAtomTypes(fitmol);
		OEMMFF94PartialCharges(fitmol);
		OEAssignBondiVdWRadii(fitmol);

		tanimoto = et.Tanimoto(fitmol);



	return tanimoto;

}
;

double NDimMolElement::PreciseValue() {
	if (x[1] == x[4] && x[2] == x[5] && x[3] == x[6]) {
		return 0;
	}

	double * newAtomsRotated =
			&(INI.getMolVariable()->getAtomsXYZ()[INI.getMolVariable()->atoms.size()
					* 3]);
	MoveAndRotate::RotateMolAccording1Axis(INI.getMolVariable()->getAtomsXYZ(),
	INI.getMolVariable()->atoms.size(), x[0], x[1], x[2], x[3], x[4], x[5],
			x[6], newAtomsRotated);

	MoveAndRotate::MolToNewPosition(newAtomsRotated,
	INI.getMolVariable()->atoms.size(), x[7], x[8], x[9]);


	////////////RESTRICCION SHAPE SIMILARITY
		double VAB = VolumeOverlap::overlapWEGA(
		INI.getMolQuery()->getAtomsXYZ(),
		INI.getMolQuery()->atoms.size(), INI.getMolQuery()->getWeightAtoms(),
		INI.getMolQuery()->getRadiusAtoms(), newAtomsRotated,
		INI.getMolVariable()->atoms.size(),
		INI.getMolVariable()->getWeightAtoms(),
		INI.getMolVariable()->getRadiusAtoms(),
		INI.getSameVanDerWaalsRadius());
		//cout << VAB << endl;
		//cout << INI.getMolQuery()->tanimoto << endl;
		//cout << INI.getMolVariable()->tanimoto << endl;
		double result2 = Tanimoto::calculateTanimotoGeneric(
		INI.getMolQuery()->tanimoto,
		INI.getMolVariable()->tanimoto, VAB);
	/////////////////////////////////////////////////

	INI.getMoleculeTemp()->updateAtomCoordinates(newAtomsRotated);

	OEGraphMol refmol, fitmol;

		oemolistream ifs;
		ifs.SetFormat(OEFormat::MOL2);
		ifs.openstring(INI.getMolQueryToString());

		OEReadMolecule(ifs, refmol);

		OEMMFFAtomTypes(refmol);
		OEMMFF94PartialCharges(refmol);
		OEAssignBondiVdWRadii(refmol);

		OEET et;
		et.SetRefMol(refmol);

		float tanimoto;
		ifs.openstring(WriteMolecule::MolToString(INI.getMoleculeTemp()));
		while (OEReadMolecule(ifs, fitmol)) {
			OEMMFFAtomTypes(fitmol);
			OEMMFF94PartialCharges(fitmol);
			OEAssignBondiVdWRadii(fitmol);

			tanimoto = et.Tanimoto(fitmol);

		}

		//return tanimoto;

		return tanimoto;
}
;

