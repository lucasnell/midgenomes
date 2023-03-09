#!/usr/bin/env python3

"""
Change the name of a file to either
  (a) have full species names if it contains an abbreviated name or
  (b) have an abbreviated species name if it contains a full name.

Notes:
- Takes as input the path to a single file.
- Will not overwrite an existing file.
- It's assumed that the abbreviation is the first character of the genus name
  (uppercase) followed by the first 5 characters of the species name (lowercase).
- It's assumed that for the full name, the genus and species are separated by
  an underscore.
- Returns an error if it doesn't find one (and only one) match from the list
  that contains both abbreviated and full species names for the chironomid
  species I'm studying.
"""

import sys
import os

if __name__ == "__main__":
    
    if (len(sys.argv) - 1) > 1:
        sys.stderr.write("Only one argument allowed. Exiting.\n")
        sys.exit(1)
    
    in_full_path = os.path.abspath(sys.argv[1])
    
    if not os.path.exists(in_full_path):
        sys.stderr.write(in_full_path + " not found\n")
        sys.exit(1)
    
    in_dir = os.path.dirname(in_full_path)
    in_file = os.path.basename(in_full_path)
    
    #' List of full species names in alphabetical order
    spp_fulls = ("Aedes_aegypti",
                 "Anopheles_stephensi",
                 "Belgica_antarctica",
                 "Chironomus_riparius",
                 "Chironomus_tentans",
                 "Clunio_marinus",
                 "Culex_quinquefasciatus",
                 "Culicoides_sonorensis",
                 "Musca_domestica",
                 "Parochlus_steinenii",
                 "Polypedilum_pembai",
                 "Polypedilum_vanderplanki",
                 "Propsilocerus_akamusi",
                 "Tanytarsus_gracilentus")
    
    #' Abbreviations for each of those names:
    spp_abbrevs = [x[0] + x.split("_")[1][:5] for x in spp_fulls]
    
    #' Match input file to abbreviated and full species names
    matches_full = [in_file.find(x) > -1 for x in spp_fulls]
    matches_abbrev = [in_file.find(x) > -1 for x in spp_abbrevs]
    
    if sum(matches_full) > 0 and sum(matches_abbrev) > 0:
        sys.stderr.write("Matches to both full names and abbreviations found.\n")
        sys.exit(1)
    elif sum(matches_full) > 0:
        if sum(matches_full) > 1:
            sys.stderr.write("More than one match to full names found.\n")
            sys.exit(1)
        #' To convert full name to abbreviation:
        str_from = [x for i, x in enumerate(spp_fulls) if matches_full[i]][0]
        str_to = [x for i, x in enumerate(spp_abbrevs) if matches_full[i]][0]
    elif sum(matches_abbrev) > 0:
        if sum(matches_abbrev) > 1:
            sys.stderr.write("More than one match to abbreviations found.\n")
            sys.exit(1)
        #' To convert abbreviation to full name:
        str_from = [x for i, x in enumerate(spp_abbrevs) if matches_abbrev[i]][0]
        str_to = [x for i, x in enumerate(spp_fulls) if matches_abbrev[i]][0]
    else:
        sys.stderr.write("No matches to full names or abbreviations found.\n")
        sys.exit(1)
    
    new_file = in_file.replace(str_from, str_to)
    new_full_path = os.path.join(in_dir, new_file)
    
    if os.path.exists(new_full_path):
        sys.stderr.write(new_full_path + " already exists.\n")
        sys.exit(1)
    print("Renaming " + in_full_path + " to " + new_full_path)
    
    os.rename(in_full_path, new_full_path)
    
    sys.exit(0)
