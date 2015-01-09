from tardisatomic.macro_atom_transition import MacroAtomTransitions

def test_create_macro_atom_structure1():
    macro_transitions = MacroAtomTransitions(None, None)
    assert len(macro_transitions.macro_atom_data.columns) == 1
    assert macro_transitions.macro_atom_data.columns[0] == 'transition_id'
    1/0
