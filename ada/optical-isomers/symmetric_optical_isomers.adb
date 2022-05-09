pragma Ada_2022;

with Text_IO;           use Text_IO;
with Ada.Command_Line;  use Ada.Command_Line;
with Ada.Environment_Variables; use Ada.Environment_Variables;

with Ada.Containers.Indefinite_Hashed_Maps;
with Ada.Strings.Hash;

procedure Symmetric_Optical_Isomers is
   
   -- This program lists all unique optical configurations of a linear
   --  molecule that has a symmetry plain in its molecular
   --  connectivity graph, like Sorbitol [1]:
   
   --         CH2OH
   --         |
   --        HC-OH     (asymmetric atom 1)
   --         |
   --      HO-CH       (asymmetric atom 2)
   --         |
   --        HC-OH     (asymmetric atom 3)
   --         |
   --        HC-OH     (asymmetric atom 4)
   --         |
   --         CH2OH
   
   -- USAGE:
   --   ./bin/symmetric_optical_isomers
   --   ./bin/symmetric_optical_isomers 2
   --   ./bin/symmetric_optical_isomers 10
   
   -- Environment variables:
   --   SYMMETRIC_OPTICAL_ISOMERS_FORMULAE
   --       when set to 1, outputs also Fisher projection formulae of
   --       the generated unique optical isomers.
   --
   --   e.g.:
   --   SYMMETRIC_OPTICAL_ISOMERS_FORMULAE=1 ./bin/symmetric_optical_isomers
   --   SYMMETRIC_OPTICAL_ISOMERS_FORMULAE=1 ./bin/symmetric_optical_isomers 2
   
   -- CAVEAT: 
   -- 
   -- This program lists *all* possible configurations and checks them
   --  for equality; since the number of configurations grows
   --  exponentially with the number of atoms, the program can be only
   --  realistically applied for short molecules:
   
   -- On an Intel(R) Xeon(R) CPU E5-1620 0 @ 3.60GHz (lscpu) with 32GB
   --  RAM, the generation of all isomers for 24 asymmetric atoms
   --  takes about 1.5 min; each next atom, the time will double, so
   --  for 30 atoms, the program will take 2**(30 - 24) = 2**6 = 64
   --  times longer, about 64 * 1.5 = 96 min., or about 1.6 hours.
   
   -- For atom counts > 30, the program with 32 bit Integer type will
   --  raise Constraint Error (for a good reason ;) ).
   
   -- Refs.:
   --
   -- [1] Wikipedia. Sorbitol. URL:
   --  https://en.wikipedia.org/wiki/Sorbitol [accessed
   --  2022-05-08T13:40+03:00, permalink:
   --  https://en.wikipedia.org/w/index.php?title=Sorbitol&oldid=1074747785]
   
   -- Each asymmetric atom can have one of the two configurations (S
   --  or R, L or D); we do not need to know which is which for just
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
      
   -- Revert molecule: list the atoms (bits, array elements) in the reverse order.
   --   
   -- This operation is equivalent to mirroring the molecule around
   --  the plane perpendicular to the molecule "axis" (the vertical
   --  axis in the above drawing):
   function Revert (Molecule : Molecule_Optical_Configuration)
                   return Molecule_Optical_Configuration is
      Reverted : Molecule_Optical_Configuration(Molecule'Range);
   begin
      for I in Molecule'Range loop
         Reverted (Molecule'Last - I + 1) := Molecule (I);
      end loop;
      return Reverted;
   end;
   
   -- Invert the molecule: change all bits to the opposite bits. This
   --  is equivalent to mirroring the molecule along the plain that is
   --  perpendicular to the Fisher projection drawing and runs along
   --  the molecule main axis (the molecule atom chain, the vertical
   --  axis in the above drawing). Combination of Revert and Invert
   --  operations is equivalent to rotation of the molecule 180 degrees
   --  around the middle axis perpendicular to the plain of the
   --  Fisher projection drawing:
   function Invert (Molecule : Molecule_Optical_Configuration)
                   return Molecule_Optical_Configuration is
      (for I in Molecule'Range => 1 - Molecule (I));
      
   -- Specify whether to print out full chemical formulae (Fisher
   --  projections) of the generated molecules:
   Print_Formulae : Boolean := False;
   
   procedure Put_Fisher_Projection
     (Molecule : in Molecule_Optical_Configuration) is
   begin
      Put_Line ("   CH2OH");
      for I in Molecule'Range loop
         Put_Line ("   |");
         if Molecule (I) = 0 then
            Put_Line ("HO-CH");
         else
            Put_Line ("  HC-OH");
         end if;
      end loop;
      Put_Line ("   |");
      Put_Line ("   CH2OH");
   end;
   
begin
   
   -- Number of atoms can be specified as a single integer number on the command line, e.g.:
   -- ./symmetric_optical_isomers 12
   if Argument_Count > 0 then
      N_Atoms := Integer'Value (Argument (1));
   end if;
   
   if Exists ("SYMMETRIC_OPTICAL_ISOMERS_FORMULAE") then
      Print_Formulae := True;
   end if;
   
   declare
      Max_Isomers : Integer := 2 ** N_Atoms;
      Molecule : Molecule_Optical_Configuration(1..N_Atoms) := (others => 0);
      Reverted : Molecule_Optical_Configuration(1..N_Atoms);
      Inverted : Molecule_Optical_Configuration(1..N_Atoms);
      
      package Atom_Configuration_Map is 
         new Ada.Containers.Indefinite_Hashed_Maps (String, Integer,
                                                    Ada.Strings.Hash, "=");
      
      use Atom_Configuration_Map;

      Observed_Molecules : Atom_Configuration_Map.Map;
   begin
      for I in 1..Max_Isomers loop
         Make_Configuration (Molecule, I);
         Inverted := Invert (Molecule);
         Reverted := Revert (Inverted);
         
         if not Contains (Observed_Molecules, To_String (Reverted)) then
            
            Insert (Observed_Molecules, To_String (Molecule), I);
            
            Put (To_String (Molecule));
            Put (" ");
            Put (To_String (Reverted));
            Put (" ");
            Put (To_String (Inverted));
            if Molecule = Reverted then
               Put (" dyad");
            end if;
            if Revert(Molecule) = Molecule then
               Put (" achiral");
            end if;
            New_Line;
            
            if Print_Formulae then
               New_Line;
               Put_Fisher_Projection (Molecule);
               New_Line;
            end if;
         end if;
      end loop;
   end;
   
end Symmetric_Optical_Isomers;
