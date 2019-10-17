def xyz(mol):
    xyz_fmt = "{0:d}\n{1:s}\n{2:s}"
    atom_line_fmt = "{0:s}{1[0]:17.6f}{1[1]:12.6f}{1[2]:12.6f}\n"

    atom_lines = ""
    for atom in mol.get_atoms():
        atom_lines += atom_line_fmt.format(atom.get_label(), atom.get_coordinate())

    nat = len(mol)
    print(xyz_fmt.format(nat, mol.get_name(), atom_lines))
