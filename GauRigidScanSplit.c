# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <stdbool.h>

/*********************************************************************
 * Assume all geomitry optimization steps converged successfullly.   *
 * This file extraxts structures of Each point in a Gaussian         *
 * rigid-scan task" output file and write to a multi-frame xyz file. *
 *********************************************************************/

int main(int argc, char const *argv[])
{
    int iarg = 0;
    char const *elements_list[] = {"", \
     "H" , "He", "Li", "Be", "B" , "C" , "N" , "O" , \
     "F" , "Ne", "Na", "Mg", "Al", "Si", "P" , "S" , \
     "Cl", "Ar", "K" , "Ca", "Sc", "Ti", "V" , "Cr", \
     "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", \
     "As", "Se", "Br", "Kr", "Rb", "Sr", "Y" , "Zr", \
     "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", \
     "In", "Sn", "Sb", "Te", "I" , "Xe", "Cs", "Ba", \
     "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", \
     "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", \
     "Ta", "W" , "Re", "Os", "Ir", "Pt", "Au", "Hg", \
     "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", \
     "Ac", "Th", "Pa", "U" , "Np", "Pu", "Am", "Cm", \
     "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", \
     "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", \
     "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};
    char ifl_name[BUFSIZ + 1] = "";
    int const num_useless_lines = 4;
    FILE *ifl = NULL;
    FILE *ofl = NULL;
    int num_frames = 0;
    int num_atoms = 0;
    int i = 0;
    int element_index = 0;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    size_t read_pos = 0;
    char buf[BUFSIZ + 1] = "";
    double ene = 0.0;
    char c = '\0';
    char coordinates_locator[] = "Standard orientation";
    char ene_type_str[11] = "";
    char *tok = NULL;

    if (argc - 1 > 1)
    {
        fprintf(stderr, "Error! Too many command arguments! At most 1 required, but got %d.\n", argc - 1);
        exit(EXIT_FAILURE);
    }
    else if (argc - 1 == 1)
    {
        strncpy(ifl_name, argv[1], BUFSIZ + 1);
    }
    else /* if (argc - 1 == 0) */
    {
        printf("Input the path of the Gaussian rigid scan output file: \n");
        fgets(ifl_name, BUFSIZ, stdin);
        if (ifl_name[strlen(ifl_name) - 1] == '\n')
        {
            ifl_name[strlen(ifl_name) - 1] = '\0';
        }
    }

    /* prepare files */
    ifl = fopen(ifl_name, "rt");
    ofl = fopen("RigidScanSplit.xyz", "rt");
    if (ofl)
    {
        fprintf(stderr, "Warning! File \"RigidScanSplit.xyz\" already exists, moving to \"RigidScanSplit.bak\".\n");
        fclose(ofl);
        ofl = fopen("RigidScanSplit.bak", "rt");
        if (ofl)
        {
            fclose(ofl);
            remove("RigidScanSplit.bak");
        }
        rename("RigidScanSplit.xyz", "RigidScanSplit.bak");
    }
    ofl = fopen("RigidScanSplit.xyz", "wt");

    /* obtain amount of atoms and coordinates locator */
    while (fgets(buf, BUFSIZ, ifl))
    {
        if (strstr(buf, "NAtoms="))
        {
            break;
        }
    }
    sscanf(strstr(buf, "NAtoms=") + strlen("NAtoms="), "%d", & num_atoms);
    while (true)
    {
        if (! fgets(buf, BUFSIZ, ifl))
        {
            fprintf(stderr, "Error! Cannot find \"Input\", \"Standard\" or \"Z-Matrix\" orientation.\n");
            fclose(ifl);
            ifl = NULL;
            fclose(ofl);
            ofl = NULL;
            exit(EXIT_FAILURE);
        }
        if (strstr(buf, coordinates_locator))
        {
            break;
        }
        if (strstr(buf, "Input orientation"))
        {
            strcpy(coordinates_locator, "Input orientation");
            break;
        }
        if (strstr(buf, "Z-Matrix orientation"))
        {
            strcpy(coordinates_locator, "Z-Matrix orientation");
            break;
        }
    }
    rewind(ifl);

    /* check energy type */
    while (true)
    {
        if (! fgets(buf, BUFSIZ, ifl))
        {
            fprintf(stderr, "Error! Cannot determine energy type.\n");
            fclose(ifl);
            ifl = NULL;
            fclose(ofl);
            ofl = NULL;
            exit(EXIT_FAILURE);
        }
        if (strstr(buf, "Energy="))
        {
            strcpy(ene_type_str, "MM");
            break;
        }
        if (strstr(buf, "SCF Done"))
        {
            break;
        }
    }
    if (strcmp(ene_type_str, "MM")) /* found "SCF Done" */
    {
        while (true)
        {
            if (! fgets(buf, BUFSIZ, ifl) || strstr(buf, "Population analysis"))
            {
                strcpy(ene_type_str, "SCF");
                break;
            }
            if (strstr(buf, "EUMP2 ="))
            {
                strcpy(ene_type_str, "MP2");
                break;
            }
            if (! strncmp(buf, " E2(", strlen(" E2(")))
            {
                strcpy(ene_type_str, "DFTPT2");
                break;
            }
            if (strstr(buf, "E(CIS/TDA)"))
            {
                strcpy(ene_type_str, "CIS/TDA");
                break;
            }
            if (strstr(buf, "E(TD-HF/TD-DFT)"))
            {
                strcpy(ene_type_str, "TD");
                break;
            }
        }
    }
    rewind(ifl);

    num_frames = 0;
    while (fgets(buf, BUFSIZ, ifl))
    {
        if (strstr(buf, coordinates_locator))
        {
            read_pos = ftell(ifl);
            /* read energy here */
            if (strstr(ene_type_str, "MM"))
            {
                while (fgets(buf, BUFSIZ, ifl))
                {
                    if (tok = strstr(buf, "Energy="))
                    {
                        break;
                    }
                }
                tok += strlen("Energy=");
                sscanf(tok, "%lg", & ene);
            }
            else if (strstr(ene_type_str, "SCF"))
            {
                while (fgets(buf, BUFSIZ, ifl))
                {
                    if (strstr(buf, "SCF Done"))
                    {
                        break;
                    }
                }
                tok = strchr(buf, '=') + strlen("=");
                sscanf(tok, "%lg", & ene);
            }
            else if (strstr(ene_type_str, "MP2"))
            {
                while (fgets(buf, BUFSIZ, ifl))
                {
                    if (tok = strstr(buf, "EUMP2"))
                    {
                        break;
                    }
                }
                tok = strchr(tok, ' ') + strlen("=");
                * strchr(tok, 'D') = 'E';
                sscanf(tok, "%lg", & ene);
            }
            else if (strstr(ene_type_str, "DFTPT2"))
            {
                while (fgets(buf, BUFSIZ, ifl))
                {
                    if (! strncmp(buf, " E2(", strlen(" E2(")))
                    {
                        break;
                    }
                }
                tok = strstr(buf, "E(");
                tok = strchr(tok, '=') + strlen("=");
                * strchr(tok, 'D') = 'E';
                sscanf(tok, "%lg", & ene);
            }
            else if (strstr(ene_type_str, "CIS/TDA"))
            {
                while (fgets(buf, BUFSIZ, ifl))
                {
                    if (tok = strstr(buf, "E(CIS/TDA)"))
                    {
                        break;
                    }
                }
                tok = strchr(buf, '=') + strlen("=");
                sscanf(tok, "%lg", & ene);
            }
            else if (strstr(ene_type_str, "TD"))
            {
                while (fgets(buf, BUFSIZ, ifl))
                {
                    if (tok = strstr(buf, "E(TD-HF/TD-DFT)"))
                    {
                        break;
                    }
                }
                tok = strchr(buf, '=') + strlen("=");
                sscanf(tok, "%lg", & ene);
            }
            else
            {
                /* should never happen */
                ;
            }
            tok = NULL;
            fseek(ifl, read_pos, SEEK_SET);
            /* read geometry here */
            ++ num_frames;
            fprintf(ofl, "%4d\n", num_atoms);
            fprintf(ofl, "frame %4d: energy = %17.10lf Hartree (\"%s\" level)\n", num_frames, ene, ene_type_str);
            for (i = 0; i < num_useless_lines; ++ i)
            {
                fgets(buf, BUFSIZ, ifl);
            }
            for (i = 0; i < num_atoms; ++ i)
            {
                fgets(buf, BUFSIZ, ifl);
                sscanf(buf, "%*d %d %*d %lg %lg %lg", & element_index, & x, & y, & z);
                fprintf(ofl, " %-2s %15s %12.8lf    %12.8lf    %12.8lf\n", elements_list[element_index], "", x, y, z);
            }
            fflush(ofl);
        }
    }

    fclose(ifl);
    ifl = NULL;
    fclose(ofl);
    ofl = NULL;

    printf("Total amount of atoms: %d\n", num_atoms);
    printf("Total amount of frames: %d\n", num_frames);
    printf("Done!\n");

    if (argc - 1 == 0)
    {
        printf("Press <ENTER> to exit ...\n");
        while((c = getchar()) != '\n' && c != EOF)
        {
            ;
        }
    }

    return 0;
}

