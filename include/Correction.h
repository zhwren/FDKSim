#ifndef CORRECTION_H
#define CORRECTION_H

class Correction
{
public:
	Correction(int iGeometry, int nDetU, int nDetV, double dUSize, double dVSize, double dFocalLength, double dSourceDetectorLength);
	~Correction();

	bool fail(void) const       {return m_fail;}
	const std::string& failMessage(void) const {return m_failMessage;}

	void init(int iGeometry, int nDetU, int nDetV, double dUSize, double dVSize, double dFocalLength, double dSourceDetectorLength);

	void ProcessSignal(double *pSino);

private:
	int m_idGeometry;

	unsigned int m_nDetU;       
	unsigned int m_nDetV; 
	double m_dDetUSize;			/*Size or Angle of detector element in array*/
	double m_dDetVSize;			/*Size of detector element in rows*/

	double m_dFocalLength;        // Focal Length, distance from source to center
	double m_dSourceDetectorLength; // Distance from source to detectors

	bool m_fail;
	std::string m_failMessage;

	double *m_dCorrectionMatix;
};



#endif
