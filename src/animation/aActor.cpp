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
	//vec3 guideGlobalT = m_Guide.getLocal2Global().m_translation;
	vec3 rootGlobalT = m_pSkeleton->getRootNode()->getGlobalTranslation();
	vec3 guideGlobalV = m_Guide.getLocal2Global() * rootGlobalT;
	// 2.	Set the y component of the guide position to 0
	m_Guide.setGlobalTranslation(vec3(guideGlobalV[0], 0.0, guideGlobalV[2]));
	// 3.	Set the global rotation of the guide joint towards the guideTarget
	vec3 target = (guideTargetPos - guideGlobalV);
	target[1] = 0.0;
	target = target.Normalize();
	m_Guide.setGlobalRotation(mat3(vec3(0.0f, 1.0f, 0.0f).Cross(target), vec3(0.f, 1.f, 0.f), target).Transpose());
	m_pSkeleton->update();
}

void AActor::solveFootIK(float leftHeight, float rightHeight, bool rotateLeft, bool rotateRight, vec3 leftNormal, vec3 rightNormal)
{
	if (!m_pSkeleton->getRootNode()) { return; }
	AJoint* leftFoot = m_pSkeleton->getJointByID(m_IKController->mLfootID);
	AJoint* rightFoot = m_pSkeleton->getJointByID(m_IKController->mRfootID);

	// TODO: 
	// The normal and the height given are in the world space

	// 1.	Update the local translation of the root based on the left height and the right height
	vec3 rootLT = m_pSkeleton->getRootNode()->getLocalTranslation();
	vec3 rootLT1 = vec3(rootLT[0], rootLT[1] + leftHeight, rootLT[2]);
	vec3 rootLT2 = vec3(rootLT[0], rootLT[1] + rightHeight, rootLT[2]);
	m_pSkeleton->getRootNode()->setLocalTranslation(Max(rootLT1, rootLT2));

	m_pSkeleton->update();

	// 2.	Update the character with Limb-based IK
	
	// Rotate Foot
	if (rotateLeft)
	{
		// Update the local orientation of the left foot based on the left normal
		ATarget t;
		vec3 lf = leftFoot->getGlobalTranslation();
		t.setGlobalTranslation(vec3(lf[0], 0.0, lf[2]));

		float angle = acos(Dot(leftNormal.Normalize(), vec3(0.0, 1.0, 0.0)));
		vec3 axis = leftNormal.Normalize().Cross(vec3(0.0, 1.0, 0.0));
		axis = (leftFoot->getLocal2Global().Inverse() * axis).Normalize();
		leftFoot->setLocalRotation(leftFoot->getLocalRotation().Rotation3D(axis, angle));

		m_IKController->IKSolver_Limb(leftFoot->getID(), t);
	}
	if (rotateRight)
	{
		ATarget t;
		vec3 rf = rightFoot->getGlobalTranslation();
		t.setGlobalTranslation(vec3(rf[0], rightHeight, rf[2]));

		float angle = acos(Dot(leftNormal.Normalize(), vec3(0.0, 1.0, 0.0)));
		vec3 axis = leftNormal.Normalize().Cross(vec3(0.0, 1.0, 0.0));
		axis = (rightFoot->getLocal2Global().Inverse() * axis).Normalize();
		rightFoot->setLocalRotation(rightFoot->getLocalRotation().Rotation3D(axis, angle));

		m_IKController->IKSolver_Limb(rightFoot->getID(), t);
	}
	m_pSkeleton->update();
}
