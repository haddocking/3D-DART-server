3D-DART DNA Modelling Server
Version 1.2

Copyright 2008-2010 Marc van Dijk
Bijvoet Center for Biomolecular Research
Utrecht University
Utrecht, The Netherlands

------------------------------------------------------------------------------------------------------------
Requirements:
This distribution has been succesfully tested on both Linux and MacOSX operating systems.
Running 3D-DART as web service requires:
- python2.5 or higher
- PHP enambles webserver (Apache has been tested succesfully).

Installation:
- Most of the variables are set automaticly 
- Specify the download location (FTP_LOCATION) in the /3DDART/server/system/Constants.py file.
- The number of served requests is logged to a count file. The location can be specified in
  in the variable SERVERCOUNTFILE in /3DDART/server/system/Contacts.py file. If you do not
  want to use this parameter that outcomment it.

------------------------------------------------------------------------------------------------------------
Changelog:
- 20/02/2010: First distribution of the service as self-contained distribution
