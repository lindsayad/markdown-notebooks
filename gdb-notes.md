# 8/4/17

Ok, currently looking in cp-demangle.c. 

There is this structure type called `d_growable_string`. It appears to have at
least two data members: `dgs->buf` and `dgs->alc`. A new memory space
(`realloc`) is created in `d_growable_string_resize` for a pointer called
`newbuf`. Then `dbs->buf` is made to point at `newbuff`. So we don't want to
free newbuf because we might want to access `dbg->buf` later.

Ok, so first commited change (cp-support.c):

- Before: 
```
==5923== LEAK SUMMARY:
==5923==    definitely lost: 734,882 bytes in 6,231 blocks
==5923==    indirectly lost: 42,581 bytes in 6 blocks
==5923==      possibly lost: 112,422 bytes in 327 blocks
==5923==    still reachable: 8,514,469 bytes in 16,790 blocks
==5923==         suppressed: 0 bytes in 0 blocks
```
- After:
```
==1748== LEAK SUMMARY:
==1748==    definitely lost: 74,226 bytes in 21 blocks
==1748==    indirectly lost: 42,581 bytes in 6 blocks
==1748==      possibly lost: 111,142 bytes in 324 blocks
==1748==    still reachable: 8,515,463 bytes in 16,791 blocks
==1748==         suppressed: 0 bytes in 0 blocks
```

Ok, we have this equality set-up through calling:
```c++
aymbol **ret = &synthsyms;
```
At the address of synthsyms must be a pointer to asymbol.

In the `elf64-x86-64.c` routine, we have this:
```c++
*ret = (asymbol *) bfd_zmalloc (size);
```
So we dereference ret once, and the result is of type `asymbol *`.

So, I know that:
```
asymbol **ret = &synthsyms;
```
and consequently, `synthsyms = *ret = (asymbol *) bfd_zmalloc (size);`.

So synthsyms is a pointer to asymbol. This is verified by the contents at the
beginning of the function `elf_read_minimal_symbols` in `elfread.c`. Synthsyms
is used in assignments to `synth_symbol_table`.

`synth_symbol_table.get()` returns a pointer to an array of pointers to
asymbols. The asymbols exist in dynamically allocated memory which must be
deallocated before the pointers are deleted.

```
      elf_symtab_read (reader, objfile, ST_SYNTHETIC, synthcount,
		       synth_symbol_table.get (), true);
```
In `elf_symtab_read`, sym is a pointer to asymbol.

Ok, implemented `xfree ((char *) synthsyms->name)` as well as `xfree
(synthsyms)` in `elfread.c`. This resulted in a new leak summary, whose xtree
can be viewed in `xtleak.kcg.30129`, shown below:

```
==30129== LEAK SUMMARY:
==30129==    definitely lost: 37,538 bytes in 15 blocks
==30129==    indirectly lost: 0 bytes in 0 blocks
==30129==      possibly lost: 111,142 bytes in 324 blocks
==30129==    still reachable: 8,512,473 bytes in 16,788 blocks
==30129==         suppressed: 0 bytes in 0 blocks
```

With two edits and a singe run of `gdb --args moltres-dbg -i
simple_diffusion.i`, the maximum resident set size is 3,520,252 kbytes. Removing
those two edits with the same run, the maximum resident set size is 6,838,984
kbytes.

12.7
17.4

The maximum resident set size is reproducible.

# 8/7/17

Now performing periodic leak checks. This is pretty cool. Breaking on
`run_command`. First leak report:
```
==16456== xtree leak report: /home/lindsayad/programming/cpp/xtleak.kcg.16456.1
==16456== LEAK SUMMARY:
==16456==    definitely lost: 37,538 (+37,538) bytes in 15 (+15) blocks
==16456==    indirectly lost: 0 (+0) bytes in 0 (+0) blocks
==16456==      possibly lost: 120,758 (+120,758) bytes in 339 (+339) blocks
==16456==    still reachable: 9,188,177 (+9,188,177) bytes in 17,269 (+17,269) blocks
==16456==         suppressed: 0 (+0) bytes in 0 (+0) blocks
==16456== Reachable blocks (those to which a pointer was found) are not shown.
==16456== To see them, add 'reachable any' args to leak_check
==16456== 
```
Note that the number is the same as my last leak check summary. Second leak
report after running my hello world program a second time through gdb:
```
==16456== xtree leak report: /home/lindsayad/programming/cpp/xtleak.kcg.16456.2
==16456== LEAK SUMMARY:
==16456==    definitely lost: 61,938 (+24,400) bytes in 22 (+7) blocks
==16456==    indirectly lost: 0 (+0) bytes in 0 (+0) blocks
==16456==      possibly lost: 120,030 (-728) bytes in 338 (-1) blocks
==16456==    still reachable: 9,188,943 (+766) bytes in 17,273 (+4) blocks
==16456==         suppressed: 0 (+0) bytes in 0 (+0) blocks
==16456== Reachable blocks (those to which a pointer was found) are not shown.
==16456== To see them, add 'reachable any' args to leak_check
==16456== 
```
Leak check number 3:
```
==16456== xtree leak report: /home/lindsayad/programming/cpp/xtleak.kcg.16456.3
==16456== LEAK SUMMARY:
==16456==    definitely lost: 86,338 (+24,400) bytes in 29 (+7) blocks
==16456==    indirectly lost: 0 (+0) bytes in 0 (+0) blocks
==16456==      possibly lost: 120,030 (+0) bytes in 338 (+0) blocks
==16456==    still reachable: 9,188,981 (+38) bytes in 17,276 (+3) blocks
==16456==         suppressed: 0 (+0) bytes in 0 (+0) blocks
==16456== Reachable blocks (those to which a pointer was found) are not shown.
==16456== To see them, add 'reachable any' args to leak_check
==16456== 
```
As I might hope, the number of definitely lost blocks increased by the same
amount.

Now, running with `gdb --args moltres-dbg -i simple_diffusion.i`. Initial leak
check before any moltres runs:
```
==16943== xtree leak report: /home/lindsayad/projects/moltres/problems/xtleak.kcg.16943.1
==16943== LEAK SUMMARY:
==16943==    definitely lost: 13,138 (+13,138) bytes in 8 (+8) blocks
==16943==    indirectly lost: 0 (+0) bytes in 0 (+0) blocks
==16943==      possibly lost: 88,887 (+88,887) bytes in 228 (+228) blocks
==16943==    still reachable: 22,994,700 (+22,994,700) bytes in 21,025 (+21,025) blocks
==16943==         suppressed: 0 (+0) bytes in 0 (+0) blocks
==16943== Reachable blocks (those to which a pointer was found) are not shown.
==16943== To see them, add 'reachable any' args to leak_check
==16943== 
```
After one moltres run (using `monitor leak_check xtleak increased` as I have for
all previous leak monitor leak checks):
```
Ignoring packet error, continuing...
==16943== xtree leak report: /home/lindsayad/projects/moltres/problems/xtleak.kcg.16943.2
==16943== LEAK SUMMARY:
==16943==    definitely lost: 452,114 (+438,976) bytes in 119 (+111) blocks
==16943==    indirectly lost: 680 (+680) bytes in 1 (+1) blocks
==16943==      possibly lost: 124,822 (+35,935) bytes in 340 (+112) blocks
==16943==    still reachable: 1,647,274,055 (+1,624,279,355) bytes in 353,328 (+332,303) blocks
==16943==         suppressed: 0 (+0) bytes in 0 (+0) blocks
==16943== Reachable blocks (those to which a pointer was found) are not shown.
==16943== To see them, add 'reachable any' args to leak_check
==16943== 
```
After two moltres runs (this time leak check conducted with `monitor leak_check
xtleak increased kinds all`):
```
Ignoring packet error, continuing...
==16943== xtree leak report: /home/lindsayad/projects/moltres/problems/xtleak.kcg.16943.3
==16943== LEAK SUMMARY:
==16943==    definitely lost: 891,090 (+438,976) bytes in 230 (+111) blocks
==16943==    indirectly lost: 1,360 (+680) bytes in 2 (+1) blocks
==16943==      possibly lost: 124,094 (-728) bytes in 339 (-1) blocks
==16943==    still reachable: 1,647,274,821 (+766) bytes in 353,332 (+4) blocks
==16943==         suppressed: 0 (+0) bytes in 0 (+0) blocks
==16943== 
```
Again the amount of definitely lost and idinrectly lost bits increases in a
reproducible way. Possibly lost actually decreases; still reachable increases by
a miniscule amount. A third run again produces the same increases in definitely
lost and indirectly lost bytes. Here's the final heap summary:
```
==16943== 
==16943== xtree memory report: /home/lindsayad/projects/moltres/problems/xtmemory.kcg.16943
==16943== HEAP SUMMARY:
==16943==     in use at exit: 1,647,963,978 bytes in 353,522 blocks
==16943==   total heap usage: 181,183,928 allocs, 180,830,406 frees, 284,463,197,787 bytes allocated
==16943== 
```


# 8/10/17

Tracking down chunkfun leak. `h` is a pointer to a struct of type `obstack`. The
chunk member of `h` is a pointer to dynamically allocated memory. 
