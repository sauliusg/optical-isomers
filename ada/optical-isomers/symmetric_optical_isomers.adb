pragma Ada_2022;

with Text_IO;           use Text_IO;
with Ada.Command_Line;  use Ada.Command_Line;

procedure Symmetric_Optical_Isomers is
   
   -- This program lists all unique optical cofgigurations of a linear
   --  molecule that has a symmetry plain in its molecular
   --  connectivity graph, like Sorbitol [1]:
   
   --         CH2OH
   --         |
   --        HC-OH     (assymetric atom 1)
   --         |
   --      HO-CH       (assymetric atom 2)
   --         |
   --        HC-OH     (assymetric atom 3)
   --         |
   --        HC-OH     (assymetric atom 4)
   --         |
   --         CH2OH
   
   -- CAVEAT: 
   -- 
   -- This program lists *all* possible configurations and checks them
   --  for equality; since the number of congugurations grows
   --  exponentially with the number of atoms, the program can be only
   --  realistically applied for very short molecules (max. 10-12
   --  atoms)... For longer molecules it will run out of memory, or
   --  overflow, or both.
   
   -- Refs.:
   --
   -- [1] Wikipedia. Sorbitol. URL:
   --  https://en.wikipedia.org/wiki/Sorbitol [accessed
   --  2022-05-08T13:40+03:00, permalink:
   --  https://en.wikipedia.org/w/index.php?title=Sorbitol&oldid=1074747785]
   
   -- Each asymmetric atom can have one of the two configurations (S
   --  or R, L or D); we do not need to know wich is whic for just
   --  listing the isomers, so we use a single bit, values '0' and
   --  '1', to specify the configuration:
   type Optical_Configuration is range 0..1;
   
   type Molecule_Optical_Configuration is array (Integer range <>) of Optical_Configuration;
   
   -- Number of asymmetric atoms to consider, with a default value:
   N_Atoms : Integer := 4;
   
   -- Make a molecular configuration array from a number; essentially
   --  fill the array with the bits from that number:
   procedure Make_Configuration (Molecule : out Molecule_Optical_Configuration; 
                                 N : Positive ) is
      I : Integer := Molecule'First;
      Bits : Natural := N - 1;
   begin
      Molecule := (others => 0);
      while Bits /= 0 loop
         Molecule (I) := Optical_Configuration (Bits mod 2);
         Bits := Bits / 2;
         I := I + 1;
      end loop;
   end;
   
   -- Convert the molecule array to a string. Useful for printing out and hashing:
   function To_String (Molecule : Molecule_Optical_Configuration) return String is
      (for I in Molecule'Range => (if Molecule (I) = 0 then '0' else '1'));
   
begin
   
   -- Number of atoms can be specified as a single integer number on the command line, e.g.:
   -- ./symmetric_optical_isomers 12
   if Argument_Count > 0 then
      N_Atoms := Integer'Value (Argument (1));
   end if;
   
   declare
      Max_Isomers : Integer := 2 ** N_Atoms;
      Molecule : Molecule_Optical_Configuration(1..N_Atoms) := (others => 0);
   begin
      for I in 1..Max_Isomers loop
         Make_Configuration (Molecule, I);
         
         Put (To_String (Molecule));
         New_Line;
      end loop;
   end;
   
end Symmetric_Optical_Isomers;
