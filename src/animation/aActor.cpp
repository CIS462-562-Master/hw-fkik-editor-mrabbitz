#include "aActor.h"

#pragma warning(disable : 4018)



/****************************************************************
*
*    	    Actor functions
*
****************************************************************/

AActor::AActor() 
{
	m_pInternalSkeleton = new ASkeleton();
	m_pSkeleton = m_pInternalSkeleton;

	m_BVHController = new BVHController();
	m_BVHController->setActor(this);

	m_IKController = new IKController();
	m_IKController->setActor(this);

	// code to update additional Actor data goes here
	resetGuide();

}

AActor::AActor(const AActor* actor)
{
	*this = *actor;
}

AActor& AActor::operator = (const AActor& actor)
{
	// Performs a deep copy
	if (&actor == this)
	{
		return *this;
	}
	m_pSkeleton = actor.m_pSkeleton;

	// code to update additional Actor data goes here


	return *this;
}

AActor::~AActor()
{
	 delete m_IKController;
	 delete m_BVHController;
	 delete m_pInternalSkeleton;

}

void AActor::clear()
{
	// looks like it is clearing more times than the number of actors.  as a result, m_pSkeleton is not defined for last case.
	m_pSkeleton->clear();  

	// code to update additional Actor data goes here
}

void AActor::update()
{
	if (!m_pSkeleton->getRootNode() )
		 return; // Nothing loaded
	else m_pSkeleton->update();

	// code to update additional Actor data goes here

}

ASkeleton* AActor::getSkeleton()
{
	return m_pSkeleton;
}

void AActor::setSkeleton(ASkeleton* pExternalSkeleton)
{
	m_pSkeleton = pExternalSkeleton;
}

void AActor::resetSkeleton()
{
	m_pSkeleton = m_pInternalSkeleton;
}

BVHController* AActor::getBVHController()
{
	return m_BVHController;
}

IKController* AActor::getIKController()
{
	return m_IKController;
}

void AActor::updateGuideJoint(vec3 guideTargetPos)
{
	if (!m_pSkeleton->getRootNode()) { return; }

	// TODO: 
	// 1.	Set the global position of the guide joint to the global position of the root joint
	// 2.	Set the y component of the guide position to 0
	// 3.	Set the global rotation of the guide joint towards the guideTarget

	// In our code, the root joint's parent is a nullptr because we assume the "world" is the parent of the root joint in previous parts.
	// While in Unity part, we can consider that the guide joint is the "world", so root->getGlobalTranslation() returns the translation w.r.t. the guide.
	// *** THE BELOW IS UNITY ***

	vec3 guidePos;
	AJoint *root = m_pSkeleton->getRootNode();
	guidePos = m_Guide.getGlobalTranslation() + m_Guide.getGlobalRotation() * root->getGlobalTranslation();
	guidePos[1] = 0.0;
	m_Guide.setGlobalTranslation(guidePos);

	guideTargetPos[1] = 0.0;
	vec3 forward_Global = vec3(0.0, 0.0, 1.0);
	vec3 forward_Guide = (guideTargetPos - m_Guide.getGlobalTranslation()).Normalize();

	double cosAlpha = Dot(forward_Guide, forward_Global) / (forward_Guide.Length() * forward_Global.Length());
	// clamping between -1 and 1
	cosAlpha = cosAlpha < -1.0 ? -1.0 : cosAlpha > 1.0 ? 1.0 : cosAlpha;
	double alpha = acos(cosAlpha);

	vec3 rotationAxis = forward_Global.Cross(forward_Guide).Normalize();

	quat q;
	q.FromAxisAngle(rotationAxis, alpha);

	m_Guide.setGlobalRotation(q.ToRotation());
}

void AActor::solveFootIK(float leftHeight, float rightHeight, bool rotateLeft, bool rotateRight, vec3 leftNormal, vec3 rightNormal)
{
	if (!m_pSkeleton->getRootNode()) { return; }
	AJoint* leftFoot = m_pSkeleton->getJointByID(m_IKController->mLfootID);
	AJoint* rightFoot = m_pSkeleton->getJointByID(m_IKController->mRfootID);

	// TODO: 
	// The normal and the height given are in the world space

	leftNormal = leftNormal.Normalize();
	rightNormal = rightNormal.Normalize();

	AJoint *root = m_pSkeleton->getRootNode();

	// 1.	Update the local translation of the root based on the left height and the right height

	vec3 root_translation = root->getLocalTranslation();
	root_translation[1] += leftHeight < rightHeight ? leftHeight : rightHeight;
	root->setLocalTranslation(root_translation);

	m_pSkeleton->update();

	// 2.	Update the character with Limb-based IK 
	ATarget target_Ground;
	vec3 position_Foot_Global;

	position_Foot_Global = leftFoot->getGlobalTranslation();
	position_Foot_Global[1] = leftHeight + 14.0;
	target_Ground.setGlobalTranslation(position_Foot_Global);
	m_IKController->IKSolver_Limb(leftFoot->getID(), target_Ground);

	position_Foot_Global = rightFoot->getGlobalTranslation();
	position_Foot_Global[1] = rightHeight + 14.0;
	target_Ground.setGlobalTranslation(position_Foot_Global);
	m_IKController->IKSolver_Limb(rightFoot->getID(), target_Ground);

	vec3 normal_Foot_Local;
	// Rotate Foot
	if (rotateLeft)
	{
		// Update the local orientation of the left foot based on the left normal

		// Global to local
		normal_Foot_Local = (leftFoot->getGlobalRotation().Transpose() * leftNormal).Normalize();

		vec3 targetVec = normal_Foot_Local;
		vec3 currentVec = vec3(0.0, 1.0, 0.0);

		double cosAlpha = Dot(targetVec, currentVec) / (targetVec.Length() * currentVec.Length());
		// clamping between -1 and 1
		cosAlpha = cosAlpha < -1.0 ? -1.0 : cosAlpha > 1.0 ? 1.0 : cosAlpha;
		double alpha = acos(cosAlpha);

		vec3 rotationAxis = currentVec.Cross(targetVec).Normalize();

		quat q;
		q.FromAxisAngle(rotationAxis, alpha);

		leftFoot->setLocalRotation(q.ToRotation());
	}
	if (rotateRight)
	{
		// Update the local orientation of the right foot based on the right normal

		// Global to local
		normal_Foot_Local = (rightFoot->getGlobalRotation().Transpose() * rightNormal).Normalize();

		vec3 targetVec = normal_Foot_Local;
		vec3 currentVec = vec3(0.0, 1.0, 0.0);

		double cosAlpha = Dot(targetVec, currentVec) / (targetVec.Length() * currentVec.Length());
		// clamping between -1 and 1
		cosAlpha = cosAlpha < -1.0 ? -1.0 : cosAlpha > 1.0 ? 1.0 : cosAlpha;
		double alpha = acos(cosAlpha);

		vec3 rotationAxis = currentVec.Cross(targetVec).Normalize();

		quat q;
		q.FromAxisAngle(rotationAxis, alpha);

		rightFoot->setLocalRotation(q.ToRotation());
	}
	m_pSkeleton->update();
}
