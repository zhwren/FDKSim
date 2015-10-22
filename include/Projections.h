#ifndef PROJECTIONS_H
#define PROJECTIONS_H

class Scanner;

class Projections
{
public:
	Projections (Scanner* pScanner, double* pData);
	Projections(Projections& rProj);
	~Projections();

	bool fail() const {return m_fail;}
	const std::string& failMessage() const {return m_failMessage;}

	int geometry() const {return m_pscanner->geometry();}
	unsigned int nDetU() const {return m_pscanner->nDetU();}
	unsigned int nDetV() const {return m_pscanner->nDetV();}
	double detIncU() const {return m_pscanner->detIncU();}
	double detIncV() const {return m_pscanner->detIncV();}

	unsigned int nView() const {return m_pscanner->nView();}

	double startView() const {return m_pscanner->startView();}
	double rotInc() const {return m_pscanner->rotInc();}
	double rotLen() const {return m_pscanner->rotLen();}

	double pitchLen() const {return m_pscanner->pitchLen();}
	double pitchInc() const {return m_pscanner->pitchInc();}

	double focalLength() const {return m_pscanner->focalLength();}
	double sourceDetectorLength() const {return m_pscanner->sourceDetectorLength();}
	double centerDetectorLength() const {return m_pscanner->centerDetectorLength();}

	double* prjData() {return m_pData;}
	Scanner* scanner() {return m_pscanner;}

private:
	bool m_fail;
	std::string m_failMessage;

	Scanner* m_pscanner;
	double* m_pData;
};


#endif
