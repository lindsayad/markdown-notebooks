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
