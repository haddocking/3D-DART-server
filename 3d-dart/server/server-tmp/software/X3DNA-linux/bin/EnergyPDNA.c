#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <fstream>
#include <iostream>
using namespace std;
char szOutForm[2048];

#include "EnergyParams.h"

// g++ -Wall -O3 EnergyPDNA.C -o EnergyPDNA.exe




//----------------------------------------------------------------------
int DoCalcsStep( const char * szFile )
{
  ifstream hOut( szFile );
  if( !hOut )
    {
      cerr << "Failed to open file `" << szFile << "`\n";
      return( EXIT_FAILURE );
    }
 
  
  char szBuffer[1024];
  float fEnergyTotal = 0.0;
  int nbBasePairs = 0;
  while( hOut.getline( szBuffer, 1024 ) )
    {
      if( strstr( szBuffer, "Number of base-pairs:" ) != NULL )
	{
	  strtok( szBuffer, ":" );
	  nbBasePairs = atoi( strtok( NULL, " \t" ) );
	}
      
      else
	
      // get the base-pair step params line in the file
      if( strstr( szBuffer, "Local base-pair step parameters" ) != NULL )
	{
	  // skip next line
	  hOut.getline( szBuffer, 1024 );

	  for( int iStep = 1; iStep < nbBasePairs; ++iStep )
	    {
	      //                 2         3         0         4         1        5
	      //     step      Shift     Slide      Rise      Tilt      Roll     Twist
	      //   1 CG/CG      0.09      0.04      3.20     -3.22      8.52     32.73
	      hOut.getline( szBuffer, 1024 );

	      /*char *szStepNo =*/ strtok( szBuffer, " \t" );
	      char *szStep   = strtok( NULL, " \t" );

	      int b1 = CEnergyParams::Base2Index( szStep[0] );
	      int b2 = CEnergyParams::Base2Index( szStep[1] );
	      int b3 = CEnergyParams::Base2Index( szStep[3] );
	      int b4 = CEnergyParams::Base2Index( szStep[4] );

	      float fTheta[6];
	      int fIndex[6] = { 2, 3, 0, 4, 1, 5 };
	      for( int iTheta = 0; iTheta < 6; ++iTheta )
		fTheta[ fIndex[ iTheta ] ] = atoi( strtok( NULL, " \t" ) );

	      float fEnergy = 0.0;
	      for( int iTheta = 0; iTheta < 6; ++iTheta )
		for( int jTheta = iTheta; jTheta < 6; ++jTheta )
		  fEnergy += CEnergyParams::m_fSpringStep[ b1 ][ b2 ][ iTheta ][ jTheta ] *
		    ( fTheta[ iTheta ] - CEnergyParams::m_fAverageStep[ b1 ][ b2 ][ iTheta ] ) *
		    ( fTheta[ jTheta ] - CEnergyParams::m_fAverageStep[ b1 ][ b2 ][ jTheta ] );

	      
	      for( int iTheta = 0; iTheta < 6; ++iTheta )
		for( int jTheta = iTheta; jTheta < 6; ++jTheta )
		  fEnergy += CEnergyParams::m_fSpringStep[ b3 ][ b4 ][ iTheta ][ jTheta ] *
		    ( fTheta[ iTheta ] - CEnergyParams::m_fAverageStep[ b3 ][ b4 ][ iTheta ] ) *
		    ( fTheta[ jTheta ] - CEnergyParams::m_fAverageStep[ b3 ][ b4 ][ jTheta ] );
	      

	      sprintf( szOutForm, "@ step %2d, %c%c/%c%c Energy = %+8.2f\n",
		       iStep, szStep[0], szStep[1], szStep[3], szStep[4], fEnergy );
	      cout << szOutForm;

	      fEnergyTotal += fEnergy;
	    }
	}
    }


  sprintf( szOutForm, "Total Step Energy = %+8.2f\n",
	   fEnergyTotal );
  cout << szOutForm;
  return 0;
}



//----------------------------------------------------------------------
int DoCalcsBase( const char * szFile )
{
  ifstream hOut( szFile );
  if( !hOut )
    {
      cerr << "Failed to open file `" << szFile << "`\n";
      return( EXIT_FAILURE );
    }
 
  
  char szBuffer[1024];
  float fEnergyTotal = 0.0;
  int nbBasePairs = 0;
  while( hOut.getline( szBuffer, 1024 ) )
    {
      if( strstr( szBuffer, "Number of base-pairs:" ) != NULL )
	{
	  strtok( szBuffer, ":" );
	  nbBasePairs = atoi( strtok( NULL, " \t" ) );
	}
      
      else
	
      // get the base-pair params line in the file
      if( strstr( szBuffer, "Local base-pair parameters" ) != NULL )
	{
	  // skip next line
	  hOut.getline( szBuffer, 1024 );

	  for( int iBase = 1; iBase <= nbBasePairs; ++iBase )
	    {
	      //              3         4         5         0         1         2
	      //  bp        Shear    Stretch   Stagger    Buckle  Propeller  Opening
	      // 1 C-G       0.28     -0.14      0.07      6.93    -17.31     -0.61

	      hOut.getline( szBuffer, 1024 );

	      /*char *szStepNo =*/ strtok( szBuffer, " \t" );
	      char *szBase   = strtok( NULL, " \t" );

	      int b = 0;
	      int b1 = CEnergyParams::Base2Index( szBase[0] );
	      int b2 = CEnergyParams::Base2Index( szBase[2] );
	      if( ( ( b1 == CEnergyParams::A ) && ( b2 == CEnergyParams::T ) ) ||
		  ( ( b1 == CEnergyParams::T ) && ( b2 == CEnergyParams::A ) ) )
		b = 0;
	      else
	      if( ( ( b1 == CEnergyParams::C ) && ( b2 == CEnergyParams::G ) ) ||
		  ( ( b1 == CEnergyParams::G ) && ( b2 == CEnergyParams::C ) ) )
		b = 1;
	      else
		{
		  cerr << "Unknown DNA base pairing: "
		       << szBase[0] << "-" << szBase[2] << "\n";
		  return( EXIT_FAILURE );
		}

	      float fTheta[6];
	      int fIndex[6] = { 3, 4, 5, 0, 1, 2 };
	      for( int iTheta = 0; iTheta < 6; ++iTheta )
		fTheta[ fIndex[ iTheta ] ] = atoi( strtok( NULL, " \t" ) );

	      float fEnergy = 0.0;
	      for( int iTheta = 0; iTheta < 6; ++iTheta )
		for( int jTheta = iTheta; jTheta < 6; ++jTheta )
		  fEnergy += CEnergyParams::m_fSpringBase[ b ][ iTheta ][ jTheta ] *
		    ( fTheta[ iTheta ] - CEnergyParams::m_fAverageBase[ b ][ iTheta ] ) *
		    ( fTheta[ jTheta ] - CEnergyParams::m_fAverageBase[ b ][ jTheta ] );
	     	      

	      sprintf( szOutForm, "@ base %2d, %c-%c Energy = %+8.2f\n",
		       iBase, szBase[0], szBase[2], fEnergy );
	      cout << szOutForm;

	      fEnergyTotal += fEnergy;
	    }
	}
    }


  sprintf( szOutForm, "Total Pair Energy = %+8.2f\n",
	   fEnergyTotal );
  cout << szOutForm;
  return 0;
}



//----------------------------------------------------------------------
int main( int argc, char * argv[] )
{
  // program credits
  cerr << "\nParts of this program are:\n\n"
       << "\tOlson WK, Gorin AA, Lu XJ, Hock LM, Zhurkin VB.\n"
       << "\t  DNA sequence-dependent deformability\n"
       << "\t  deduced from protein-DNA crystal complexes.\n"
       << "\t  Proc Natl Acad Sci U S A. 1998 Sep 15;95(19):11163-8.\n\n"
       << "\tLankas F, Sponer J, Langowski J, Cheatham TE 3rd.\n"
       << "\t  DNA deformability at the base pair level.\n"
       << "\t  J Am Chem Soc. 2004 Apr 7;126(13):4124-5.\n\n";


  // check arguments
  if( argc < 2 )
    {
      cerr << "Usage:\n\t" << argv[0] << " [switches] 3DNA.out\n"
	   << "\t[switches] are:\n"
	   << "\t -s for step energy calcs (default).\n"
	   << "\t -b for base energy calcs.\n";
      return( EXIT_FAILURE );
    }


  int iOffset = 1;
  int (*pfnDoCalcs)( const char * szFile );
  pfnDoCalcs = DoCalcsStep;  // default calcs


  // parse arguments
  for( ; iOffset < argc; ++iOffset )
    {
      if( argv[ iOffset ][ 0 ] != '-' )
	break;

      switch( argv[ iOffset ][ 1 ] )
	{
	case 'b':
	case 'B':
	  pfnDoCalcs = DoCalcsBase;
	  break;
	  
	case 's':
	case 'S':
	default:
	  pfnDoCalcs = DoCalcsStep;
	  break;
	}
    }


  return pfnDoCalcs( argv[ iOffset ] );
}

