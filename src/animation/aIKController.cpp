#include "aIKController.h"
#include "aActor.h"
#include <Eigen\Dense>

#pragma warning (disable : 4018)

int IKController::gIKmaxIterations = 5;
double IKController::gIKEpsilon = 0.1;

// AIKchain class functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////
AIKchain::AIKchain()
{
	mWeight0 = 0.1;
}

AIKchain::~AIKchain()
{

}

AJoint* AIKchain::getJoint(int index) 
{ 
	return mChain[index]; 
}

void AIKchain::setJoint(int index, AJoint* pJoint) 
{ 
	mChain[index] = pJoint; 
}

double AIKchain::getWeight(int index) 
{ 
	return mWeights[index]; 
}

void AIKchain::setWeight(int index, double weight) 
{ 
	mWeights[index] = weight; 
}

int AIKchain::getSize() 
{ 
	return mChain.size(); 
}

std::vector<AJoint*>& AIKchain::getChain() 
{ 
	return mChain; 
}

std::vector<double>& AIKchain::getWeights() 
{ 
	return mWeights; 
}

void AIKchain::setChain(std::vector<AJoint*> chain) 
{
	mChain = chain; 
}

void AIKchain::setWeights(std::vector<double> weights) 
{ 
	mWeights = weights; 
}

// AIKController class functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////

IKController::IKController()
{
	m_pActor = NULL;
	m_pSkeleton = NULL;
	mvalidLimbIKchains = false;
	mvalidCCDIKchains = false;

	// Limb IK
	m_pEndJoint = NULL;
	m_pMiddleJoint = NULL;
	m_pBaseJoint = NULL;
	m_rotationAxis = vec3(0.0, 1.0, 0.0);

	ATransform desiredTarget = ATransform();
	mTarget0.setLocal2Parent(desiredTarget);  // target associated with end joint
	mTarget1.setLocal2Parent(desiredTarget);  // optional target associated with middle joint - used to specify rotation of middle joint about end/base axis
	mTarget0.setLocal2Global(desiredTarget);
	mTarget1.setLocal2Global(desiredTarget);

	//CCD IK
	mWeight0 = 0.1;  // default joint rotation weight value

}

IKController::~IKController()
{
}

ASkeleton* IKController::getSkeleton()
{
	return m_pSkeleton;
}

const ASkeleton* IKController::getSkeleton() const
{
	return m_pSkeleton;
}

ASkeleton* IKController::getIKSkeleton()
{
	return &mIKSkeleton;
}

const ASkeleton* IKController::getIKSkeleton() const
{
	return &mIKSkeleton;
}

AActor* IKController::getActor()
{
	return m_pActor;
}

void IKController::setActor(AActor* actor)

{
	m_pActor = actor;
	m_pSkeleton = m_pActor->getSkeleton();
}


AIKchain IKController::createIKchain(int endJointID, int desiredChainSize, ASkeleton* pSkeleton)
{
	// TODO: given the end joint ID and the desired length of the IK chain, 
	// add the corresponding skeleton joint pointers to the AIKChain "chain" data member starting with the end joint
	// desiredChainSize = -1 should create an IK chain of maximum length (where the last chain joint is the joint before the root joint)
	// also add weight values to the associated AIKChain "weights" data member which can be used in a CCD IK implemention

	if (desiredChainSize == -1)
	{
		desiredChainSize = INT_MAX;
	}

	AIKchain ikChain;

	if (desiredChainSize > 0)
	{
		std::vector<AJoint*> chain;
		std::vector<double> weights;

		int counter = 0;
		AJoint *rootNode = pSkeleton->getRootNode();
		AJoint *joint = pSkeleton->getJointByID(endJointID);

		while (counter < desiredChainSize)
		{
			chain.push_back(joint);
			weights.push_back(0.1);
			++counter;

			joint = joint->getParent();
			if (joint == rootNode)
			{
				break;
			}
		}
		ikChain.setChain(chain);
		ikChain.setWeights(weights);
	}

	return ikChain;
}



bool IKController::IKSolver_Limb(int endJointID, const ATarget& target)
{
	// Implements the analytic/geometric IK method assuming a three joint limb  

	// copy transforms from base skeleton
	mIKSkeleton.copyTransforms(m_pSkeleton);

	if (!mvalidLimbIKchains)
	{
		mvalidLimbIKchains = createLimbIKchains();
		if (!mvalidLimbIKchains) { return false; }
	}

	vec3 desiredRootPosition;

	// The joint IDs are not const, so we cannot use switch here
	if (endJointID == mLhandID)
	{
		mLhandTarget = target;
		computeLimbIK(mLhandTarget, mLhandIKchain, -axisY, &mIKSkeleton);
	}
	else if (endJointID == mRhandID)
	{
		mRhandTarget = target;
		computeLimbIK(mRhandTarget, mRhandIKchain, axisY, &mIKSkeleton);
	}
	else if (endJointID == mLfootID)
	{
		mLfootTarget = target;
		computeLimbIK(mLfootTarget, mLfootIKchain, axisX, &mIKSkeleton);
	}
	else if (endJointID == mRfootID)
	{
		mRfootTarget = target;
		computeLimbIK(mRfootTarget, mRfootIKchain, axisX, &mIKSkeleton);
	}
	else if (endJointID == mRootID)
	{
		desiredRootPosition = target.getGlobalTranslation();
		mIKSkeleton.getJointByID(mRootID)->setLocalTranslation(desiredRootPosition);
		mIKSkeleton.update();
		computeLimbIK(mLhandTarget, mLhandIKchain, -axisY, &mIKSkeleton);
		computeLimbIK(mRhandTarget, mRhandIKchain, axisY, &mIKSkeleton);
		computeLimbIK(mLfootTarget, mLfootIKchain, axisX, &mIKSkeleton);
		computeLimbIK(mRfootTarget, mRfootIKchain, axisX, &mIKSkeleton);
	}
	else
	{
		mIKchain = createIKchain(endJointID, 3, &mIKSkeleton);
		computeLimbIK(target, mIKchain, axisY, &mIKSkeleton);
	}

	// update IK Skeleton transforms
	mIKSkeleton.update();

	// copy IK skeleton transforms to main skeleton
	m_pSkeleton->copyTransforms(&mIKSkeleton);

	return true;
}



int IKController::createLimbIKchains()
{
	bool validChains = false;
	int desiredChainSize = 3;

	// create IK chains for Lhand, Rhand, Lfoot and Rfoot 
	mLhandIKchain = createIKchain(mLhandID, desiredChainSize, &mIKSkeleton);
	mRhandIKchain = createIKchain(mRhandID, desiredChainSize, &mIKSkeleton);
	mLfootIKchain = createIKchain(mLfootID, desiredChainSize, &mIKSkeleton);
	mRfootIKchain = createIKchain(mRfootID, desiredChainSize, &mIKSkeleton);
	
	if (mLhandIKchain.getSize() == 3 && mRhandIKchain.getSize() == 3 && mLfootIKchain.getSize() == 3 && mRfootIKchain.getSize() == 3)
	{
		validChains = true;
		
		// initalize end joint target transforms for Lhand, Rhand, Lfoot and Rfoot based on current position and orientation of joints
		mIKSkeleton.copyTransforms(m_pSkeleton);
		mLhandTarget.setLocal2Global(mIKSkeleton.getJointByID(mLhandID)->getLocal2Global());
		mRhandTarget.setLocal2Global(mIKSkeleton.getJointByID(mRhandID)->getLocal2Global());
		mLfootTarget.setLocal2Global(mIKSkeleton.getJointByID(mLfootID)->getLocal2Global());
		mRfootTarget.setLocal2Global(mIKSkeleton.getJointByID(mRfootID)->getLocal2Global());
	}

	return validChains;
}




int IKController::computeLimbIK(ATarget target, AIKchain& IKchain, const vec3 midJointAxis, ASkeleton* pIKSkeleton)
{
	// TODO: Implement the analytic/geometric IK method assuming a three joint limb  
	// The actual position of the end joint should match the target position within some episilon error 
	// the variable "midJointAxis" contains the rotation axis for the middle joint

	//// Limb (Geometric) IK variables
	//ATarget mTarget0;       // target associated with end joint
	//ATarget mTarget1;       // target associated with middle joint - used to specify rotation of middle joint about end/base axis  
	//AJoint* m_pEndJoint;    // Joint 2
	//AJoint* m_pMiddleJoint; // Joint 1
	//AJoint* m_pBaseJoint;   // Joint 0
	//double m_length01;
	//double m_length12;
	//vec3 m_rotationAxis;    // axis of rotation for middle joint.  Values can be axisX, axisY or axisZ

	m_pEndJoint = IKchain.getJoint(0);
	m_pMiddleJoint = IKchain.getJoint(1);
	m_pBaseJoint = IKchain.getJoint(2);

	// Translation of joint w.r.t. it's parent => getLocalTranslation()
	m_length01 = m_pMiddleJoint->getLocalTranslation().Length();
	m_length12 = m_pEndJoint->getLocalTranslation().Length();
	double Rd = (target.getGlobalTranslation() - m_pBaseJoint->getGlobalTranslation()).Length();

	double cosPhi = (m_length01 * m_length01 + m_length12 * m_length12 - Rd * Rd) / (2.0 * m_length01 * m_length12);
	// clamping between -1 and 1
	cosPhi = cosPhi < -1.0 ? -1.0 : cosPhi > 1.0 ? 1.0 : cosPhi;

	double theta2 = M_PI - acos(cosPhi);

	m_pMiddleJoint->setLocalRotation(mat3::Rotation3D(midJointAxis, theta2));
	m_pMiddleJoint->updateTransform();

	vec3 targetVec = target.getGlobalTranslation() - m_pBaseJoint->getGlobalTranslation();
	vec3 currentVec = m_pEndJoint->getGlobalTranslation() - m_pBaseJoint->getGlobalTranslation();

	double cosAlpha = Dot(targetVec, currentVec) / (targetVec.Length() * currentVec.Length());
	// clamping between -1 and 1
	cosAlpha = cosAlpha < -1.0 ? -1.0 : cosAlpha > 1.0 ? 1.0 : cosAlpha;
	double alpha = acos(cosAlpha);

	vec3 rotationAxis = currentVec.Cross(targetVec).Normalize();

	quat q;
	q.FromAxisAngle(rotationAxis, alpha);

	m_pBaseJoint->setLocalRotation(m_pBaseJoint->getLocalRotation() * q.ToRotation());
	m_pBaseJoint->updateTransform();

	return true;
}

bool IKController::IKSolver_CCD(int endJointID, const ATarget& target)
{
	// Implements the CCD IK method assuming a three joint limb 

	bool validChains = false;

	if (!mvalidCCDIKchains)
	{
		mvalidCCDIKchains = createCCDIKchains();
		assert(mvalidCCDIKchains);
	}

	// copy transforms from base skeleton
	mIKSkeleton.copyTransforms(m_pSkeleton);

	vec3 desiredRootPosition;


	if (endJointID == mLhandID)
	{
		mLhandTarget = target;
		computeCCDIK(mLhandTarget, mLhandIKchain, &mIKSkeleton);
	}
	else if (endJointID == mRhandID)
	{
		mRhandTarget = target;
		computeCCDIK(mRhandTarget, mRhandIKchain, &mIKSkeleton);
	}
	else if (endJointID == mLfootID)
	{
		mLfootTarget = target;
		computeCCDIK(mLfootTarget, mLfootIKchain, &mIKSkeleton);
	}
	else if (endJointID == mRfootID)
	{
		mRfootTarget = target;
		computeCCDIK(mRfootTarget, mRfootIKchain, &mIKSkeleton);
	}
	else if (endJointID == mRootID)
	{
		desiredRootPosition = target.getGlobalTranslation();
		mIKSkeleton.getJointByID(mRootID)->setLocalTranslation(desiredRootPosition);
		mIKSkeleton.update();
		computeCCDIK(mLhandTarget, mLhandIKchain, &mIKSkeleton);
		computeCCDIK(mRhandTarget, mRhandIKchain, &mIKSkeleton);
		computeCCDIK(mLfootTarget, mLfootIKchain, &mIKSkeleton);
		computeCCDIK(mRfootTarget, mRfootIKchain, &mIKSkeleton);
	}
	else
	{
		mIKchain = createIKchain(endJointID, -1, &mIKSkeleton);
		computeCCDIK(target, mIKchain, &mIKSkeleton);
	}

	// update IK Skeleton transforms
	mIKSkeleton.update();

	// copy IK skeleton transforms to main skeleton
	m_pSkeleton->copyTransforms(&mIKSkeleton);

	return true;
}

int IKController::createCCDIKchains()
{
	bool validChains = false;

	int desiredChainSize = -1;  // default of -1 creates IK chain of maximum length from end joint to child joint of root


	// create IK chains for Lhand, Rhand, Lfoot and Rfoot 
	mLhandIKchain = createIKchain(mLhandID, desiredChainSize, &mIKSkeleton);
	mRhandIKchain = createIKchain(mRhandID, desiredChainSize, &mIKSkeleton);
	mLfootIKchain = createIKchain(mLfootID, desiredChainSize, &mIKSkeleton);
	mRfootIKchain = createIKchain(mRfootID, desiredChainSize, &mIKSkeleton);

	if (mLhandIKchain.getSize() > 1 && mRhandIKchain.getSize() > 1 && mLfootIKchain.getSize() > 1 && mRfootIKchain.getSize() > 1)
	{
		validChains = true;

		// initalize end joint target transforms for Lhand, Rhand, Lfoot and Rfoot based on current position and orientation of joints
		mIKSkeleton.copyTransforms(m_pSkeleton);
		mLhandTarget.setLocal2Global(mIKSkeleton.getJointByID(mLhandID)->getLocal2Global());
		mRhandTarget.setLocal2Global(mIKSkeleton.getJointByID(mRhandID)->getLocal2Global());
		mLfootTarget.setLocal2Global(mIKSkeleton.getJointByID(mLfootID)->getLocal2Global());
		mRfootTarget.setLocal2Global(mIKSkeleton.getJointByID(mRfootID)->getLocal2Global());
	}

	return validChains;
}


int IKController::computeCCDIK(ATarget target, AIKchain& IKchain, ASkeleton* pIKSkeleton)
{

	// TODO: Implement CCD IK  
	// The actual position of the end joint should match the desiredEndPos within some episilon error 


	//Hint:
	// 1. compute axis and angle for a joint in the IK chain (distal to proximal) in global coordinates
	// 2. once you have the desired axis and angle, convert axis to local joint coords 
	// 3. compute desired change to local rotation matrix
	// 4. set local rotation matrix to new value
	// 5. update transforms for joint and all children

	AJoint *endJoint = IKchain.getJoint(0);
	AJoint *currJoint;
	vec3 errorVec;

	for (int i = 0; i < gIKmaxIterations; ++i)
	{
		for (int j = 1; j < IKchain.getSize(); ++j)
		{
			errorVec = target.getGlobalTranslation() - endJoint->getGlobalTranslation();
			if (errorVec.Length() < gIKEpsilon)
			{
				return true;
			}

			currJoint = IKchain.getJoint(j);
			vec3 r_jn = endJoint->getGlobalTranslation() - currJoint->getGlobalTranslation();

			double angle = IKchain.getWeight(j) * r_jn.Cross(errorVec).Length() / (Dot(r_jn, r_jn) + Dot(r_jn, errorVec));
			vec3 axis_Global = r_jn.Cross(errorVec).Normalize();
			vec3 axis_Local = (currJoint->getGlobalRotation().Transpose() * axis_Global).Normalize();

			currJoint->setLocalRotation(currJoint->getLocalRotation() * mat3::Rotation3D(axis_Local, angle));
			currJoint->updateTransform();
		}
	}

	return true;
}


bool IKController::IKSolver_PseudoInv(int endJointID, const ATarget& target)
{
	// TODO: Implement Pseudo Inverse-based IK  
	// The actual position of the end joint should match the target position after the skeleton is updated with the new joint angles

	if (!mvalidPseudoIKchains)
	{
		mvalidPseudoIKchains = createPseudoIKchains();
		if (!mvalidPseudoIKchains) { return false; }
	}

	// copy transforms from base skeleton
	mIKSkeleton.copyTransforms(m_pSkeleton);

	vec3 desiredRootPosition;

	if (endJointID == mLhandID)
	{
		mLhandTarget = target;
		computePseudoIK(mLhandTarget, mLhandIKchain, &mIKSkeleton);
	}
	else if (endJointID == mRhandID)
	{
		mRhandTarget = target;
		computePseudoIK(mRhandTarget, mRhandIKchain, &mIKSkeleton);
	}
	else if (endJointID == mLfootID)
	{
		mLfootTarget = target;
		computePseudoIK(mLfootTarget, mLfootIKchain, &mIKSkeleton);
	}
	else if (endJointID == mRfootID)
	{
		mRfootTarget = target;
		computePseudoIK(mRfootTarget, mRfootIKchain, &mIKSkeleton);
	}
	else if (endJointID == mRootID)
	{
		desiredRootPosition = target.getGlobalTranslation();
		mIKSkeleton.getJointByID(mRootID)->setLocalTranslation(desiredRootPosition);
		mIKSkeleton.update();
		computePseudoIK(mLhandTarget, mLhandIKchain, &mIKSkeleton);
		computePseudoIK(mRhandTarget, mRhandIKchain, &mIKSkeleton);
		computePseudoIK(mLfootTarget, mLfootIKchain, &mIKSkeleton);
		computePseudoIK(mRfootTarget, mRfootIKchain, &mIKSkeleton);
	}
	else
	{
		mIKchain = createIKchain(endJointID, 3, &mIKSkeleton);
		//mIKchain = createIKchain(endJointID, -1, &mIKSkeleton);
		computePseudoIK(target, mIKchain, &mIKSkeleton);
	}

	// update IK Skeleton transforms
	mIKSkeleton.update();

	// copy IK skeleton transforms to main skeleton
	m_pSkeleton->copyTransforms(&mIKSkeleton);

	return true;
}

int IKController::createPseudoIKchains()
{
	bool validChains = false;
	int desiredChainSize = 3;
	//int desiredChainSize = -1;  // default of -1 creates IK chain of maximum length from end joint to child joint of root

	// create IK chains for Lhand, Rhand, Lfoot and Rfoot 
	mLhandIKchain = createIKchain(mLhandID, desiredChainSize, &mIKSkeleton);
	mRhandIKchain = createIKchain(mRhandID, desiredChainSize, &mIKSkeleton);
	mLfootIKchain = createIKchain(mLfootID, desiredChainSize, &mIKSkeleton);
	mRfootIKchain = createIKchain(mRfootID, desiredChainSize, &mIKSkeleton);

	if (mLhandIKchain.getSize() == 3 && mRhandIKchain.getSize() == 3 && mLfootIKchain.getSize() == 3 && mRfootIKchain.getSize() == 3)
	//if (mLhandIKchain.getSize() > 1 && mRhandIKchain.getSize() > 1 && mLfootIKchain.getSize() > 1 && mRfootIKchain.getSize() > 1)
	{
		validChains = true;

		// initalize end joint target transforms for Lhand, Rhand, Lfoot and Rfoot based on current position and orientation of joints
		mIKSkeleton.copyTransforms(m_pSkeleton);
		mLhandTarget.setLocal2Global(mIKSkeleton.getJointByID(mLhandID)->getLocal2Global());
		mRhandTarget.setLocal2Global(mIKSkeleton.getJointByID(mRhandID)->getLocal2Global());
		mLfootTarget.setLocal2Global(mIKSkeleton.getJointByID(mLfootID)->getLocal2Global());
		mRfootTarget.setLocal2Global(mIKSkeleton.getJointByID(mRfootID)->getLocal2Global());
	}

	return validChains;
}

int IKController::computePseudoIK(ATarget target, AIKchain& IKchain, ASkeleton* pIKSkeleton)
{
	AJoint *endJoint = IKchain.getJoint(0);
	AJoint *currJoint;
	// Need to figure this out next
	Eigen::MatrixXd L(3, 3);

	Eigen::MatrixXd B(3, 3);
	vec3 r_0_jn;
	mat3 R_0_j;
	vec3 b_x, b_y, b_z;

	std::vector<Eigen::MatrixXd> J_Vec;

	double cosX, cosY, sinX, sinY;

	for (int i = 1; i < IKchain.getSize(); ++i)
	{
		currJoint = IKchain.getJoint(i);
		vec3 eulerAngles;
		currJoint->getLocalRotation().ToEulerAngles(mat3::ZYX, eulerAngles);

		std::cout << eulerAngles[0] << std::endl << eulerAngles[1] << std::endl << eulerAngles[2] << std::endl;

		cosX = cos(eulerAngles[0]);
		cosY = cos(eulerAngles[1]);
		sinX = sin(eulerAngles[0]);
		sinY = sin(eulerAngles[1]);

		L <<	1.0,  0.0	, -sinY,
				0.0,  cosX	,  sinX * cosY,
				0.0, -sinX	,  cosX * cosY;

		std::cout << L << std::endl;

		r_0_jn = endJoint->getGlobalTranslation() - currJoint->getGlobalTranslation();
		R_0_j = currJoint->getGlobalRotation().Transpose();

		std::cout << R_0_j << std::endl;
		std::cout << r_0_jn << std::endl;

		b_x = R_0_j[0].Cross(r_0_jn);
		b_y = R_0_j[1].Cross(r_0_jn);
		b_z = R_0_j[2].Cross(r_0_jn);

		B <<	b_x[0], b_y[0], b_z[0],
				b_x[1], b_y[1], b_z[1],
				b_x[2], b_y[2], b_z[2];

		std::cout << B << std::endl;

		J_Vec.push_back(B * L);
	}

	Eigen::MatrixXd J(3, 0);
	Eigen::MatrixXd currentJ;

	for (int i = J_Vec.size() - 1; i >= 0; --i)
	{
		currentJ = J_Vec[i];
		J.conservativeResize(J.rows(), J.cols() + currentJ.cols());
		J.rightCols(currentJ.cols()) = currentJ;

		std::cout << currentJ << std::endl;
		std::cout << J << std::endl;
	}

	// Left Pseudo Inverse
	// We are doing limb joints, so we only have 2 joints to rotate
	// if (J_Vec.size() < 4)
	Eigen::MatrixXd JPseduo;// = (J.transpose() * J).inverse() * J.transpose();	// should be 2 X 3
	Eigen::MatrixXd J0, J1;
	J0 = J.transpose() * J;
	std::cout << J0 << std::endl << std::endl;
	J1 = J0.inverse();
	std::cout << J1 << std::endl << std::endl;
	JPseduo = J1 * J.transpose();
	std::cout << JPseduo << std::endl << std::endl;

	vec3 endJoint_Translation = target.getGlobalTranslation() - endJoint->getGlobalTranslation();

	Eigen::Vector3d v;
	v << endJoint_Translation[0], endJoint_Translation[1], endJoint_Translation[2];
	std::cout << v << std::endl << std::endl;

	vec3 eulerAnglesShoulder;
	IKchain.getJoint(2)->getLocalRotation().ToEulerAngles(mat3::ZYX, eulerAnglesShoulder);
	vec3 eulerAnglesElbow;
	IKchain.getJoint(1)->getLocalRotation().ToEulerAngles(mat3::ZYX, eulerAnglesElbow);

	Eigen::VectorXd angles_Previous(6);
	angles_Previous <<	eulerAnglesShoulder[0], eulerAnglesShoulder[1], eulerAnglesShoulder[2],
						eulerAnglesElbow[0],	eulerAnglesElbow[1],	eulerAnglesElbow[2];

	std::cout << angles_Previous << std::endl << std::endl;

	Eigen::VectorXd result = angles_Previous + (JPseduo * v);

	eulerAnglesShoulder = vec3(result(0), result(1), result(2));
	eulerAnglesElbow = vec3(result(3), result(4), result(5));

	return true;
}

bool IKController::IKSolver_Other(int endJointID, const ATarget& target)
{
	// TODO: Put Optional IK implementation or enhancements here
	 
	return true;
}