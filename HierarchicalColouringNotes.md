**Hierarchical Colouring - notes**

**The Hierarchical Colouring approach**

Hierarchical colouring assigns a float value ("HiCol") to each tile on the basis of its hierarchy sequence (H,T,P,F). HiCol is used to assign a colour to the tile therefore attempting to highlight the global-local structure of tiles and metatiles.

It works recursively: each metatile inherits a HiCol from its parent and transmits a HiCol to its children, calculated from the inherited colour and its own colour.

The colour propagation recursion is implemented during the drawing construction, in which each metatile generates its children metatiles: 

- **Recursion** **setup**: A value "HiCol is assigned to each Metatile type (H,T,P,F) in range [-1,1], the values vector is named Motif. Similarly a values is assigned to each hat by clusters (H1,H,T,P,F), used at level 0,  e.g. { H1: -0.8, H: 0.8, T: -0.0, P: -0.8, F: 0.8 }.
- **Recursion** **initialization**: the starting inherited HiCol is 0; and HiCol=0 is assigned to the initial top level metatile in order to avoid an overall image bias.
- **Recursion**: during the recursive drawing process each metatile tile passes its HiCol plus the inherited one from its hierarchy to its children.
- **Recursion** **end**: the recursion ends at Level 0 (hats), where the drawPolygon() function is performed. The resulting HiCol value is normalized to the range [-1,1] by considering the number of levels (HiContrast() function), and mapped to an RGB colour through a Palette (ColMap() function).



**Implemented features**

- Selection of Motif vectors in menu
- Selection of Palettes in menu, sequential and divergent types
- "Save png" file has been modified to produce filenames indicating the main used parameters, and therefore compare the many possible variations.

You can test your own motiv or palette by adding them in MotifList and PaletteList (the following logic is automatic). Please let us know if you find any interesting combination…!

**Enhancement ideas**

- Two variables are defined to create a different motif at level 0 / hats (HatMotif) vs other levels / metatiles (MetaMotif), in the current implementation HatMotif=MetaMotif. A separated motif choice could be activated by adding a menu.
- HiContrast()  manages a parameter (prgContrast) which changes the contrast distribution among levels: range (-1: highlights top metatiles, 1: highlights hats), 0: neutral. It could be added to the menu.

**Contribution**

2023-04 Hierarchical colouring approach created by Santo Leonardo from hatviz SW <https://github.com/isohedral/hatviz> (Copyright (c) 2023, Craig S. Kaplan).

Mention Santo Leonardo ONLY referring to the Hierarchical Colouring approach.


