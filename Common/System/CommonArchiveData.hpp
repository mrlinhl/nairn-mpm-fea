/********************************************************************************
    CommonArchiveData.hpp
    nairn-mpm-fea
    
    Created by John Nairn on Jan 15, 2006
    Copyright (c) 2006 John A. Nairn, All rights reserved.
	
	Dependencies
		none
********************************************************************************/

#ifndef _COMMONARCHIVEDATA_

#define _COMMONARCHIVEDATA_

class CommonArchiveData
{
    public:
		static char folderDelim;
	
        //  Constructors and Destructor
        CommonArchiveData();
		virtual ~CommonArchiveData();
		
		// abstract methods
	
		// methods
		void ArchiveNodalPoints(int);
		void ArchiveElements(int);

		// accessors
		void DefineArchiveRoot(const char *);
		bool SetArchiveRoot(bool);
		char *GetArchiveRoot(void);
		void SetArchiveMesh(bool);
		void SetInputDirPath(const char *,bool);
		char *ExpandOutputPath(const char *);
		void GetFilePath(char *,const char *);
		void GetFilePathNum(char *,const char *,int);

	protected:
		char *inputDir,*archiveRoot,*globalFile,*archiveParent,*archiveDosRoot;
		char *outputDir;
		bool hasSetArchiveRoot;
		bool archiveMesh,forceUnique;
};

#endif
