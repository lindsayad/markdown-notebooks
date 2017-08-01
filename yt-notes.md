# 7/28/17

Working on issue \#1511. So I've run into one issue so far. When
`_read_fluid_fields` is called on line 1323 of `data_containers.py`, it
incorrectly return `('all', 'u')` as a `read_fluid` as opposed to a field that
needs to be generated.

Ok, `Dataset.field_dependencies` is an object I need to investigate next. How
and when does it get filled? How and when does it get updated when a derived
field is added?
