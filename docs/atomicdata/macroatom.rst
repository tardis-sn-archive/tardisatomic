.. _macro_atom_description:

***************
Macro Atom Data
***************

.. currentmodule:: tardisatomic.macro_atom_transition

The macro atom is described in detail in :cite:`2002A&A...384..725L`. The basic principal is that when an energy packet
is absorbed that the macro atom is on a certain level. A set of probabilities govern the next steps. For example, the P\ :sub:`up`,
P\ :sub:`down` and P\ :sub:`down emission` being the probability for going to a higher level, a lower level and a lower
level and emitting a photon while doing this respectively (see Figure 1 in :cite:`2002A&A...384..725L` ).

The macro atom is the most complex idea to implement as a data structure. The setup is done here in `~tardisatomic` for the coefficients
that are known beforehand.

For each level we look at the line list to see what transitions are possible
for each transition class (:ref:`transition_classes`). For each of the transition
probability classes we generate one `~pandas.DataFrame`. Each level contains a set of probabilities, The first
part of each set contains the upwards probabilities (internal upward), the second part the downwards probabilities
(internal downward), and the last part is the downward and emission probability.


The second array is for book-keeping it has exactly the length as levels (with an example for the Si II level 15):

+--------+------------------+------------+----------------+-----------------+
|Level ID| Probability index|N\ :sub:`up`| N\ :sub:`down` | N\ :sub:`total` |
+========+==================+============+================+=================+
|14001015| ???              |17          | 5              | 17 + 5*2 = 27   |
+--------+------------------+------------+----------------+-----------------+




.. _transition_classes:

Transition Categories
=====================

There are several transition categories characterised by different identification numbers:


+----------------------------+-------------+--------------------------+------------------------------+
|Name                        |Transition ID|Class                     | Section                      |
+----------------------------+-------------+--------------------------+------------------------------+
|Emission Down               |1            |:class:`PEmissionDown`    |:ref:`rad_int_bound_bound`    |
+----------------------------+-------------+--------------------------+------------------------------+

.. _rad_int_bound_bound:

Radiative and Internal Bound-Bound
----------------------------------

The three internal and radiative bound-bound transition classes P\ :sub:`up` (``id=1``),
P\ :sub:`down` (``id=2``) and P\ :sub:`down emission` (``id=3``; :class:`PEmissionDown`) are calculated

We now will calculate the transition probabilites, using the radiative rates
in Equation 20, 21, and 22 in :cite:`2002A&A...384..725L`. Then we calculate the
downward emission probability from Equation 5, the downward and upward internal
transition probabilities in :cite:`2003A&A...403..261L`.

.. math::
    p_\textrm{emission down}&= {\cal R}_{\textrm{i}\rightarrow\textrm{lower}}\,(\epsilon_\textrm{upper} - \epsilon_\textrm{lower}) / {\cal D}_{i}\\
    p_\textrm{internal down}&= {\cal R}_{\textrm{i}\rightarrow\textrm{lower}}\,\epsilon_\textrm{lower}/{\cal D}_{i}\\,
    p_\textrm{internal up}&={\cal R}_{\textrm{i}\rightarrow\textrm{upper}}\,\epsilon_{i}/{\cal D}_{i}\\,

where :math:`i` is the current level, :math:`\epsilon` is the energy of the level,
and :math:`{\cal R}` is the radiative rates.


We ignore the probability to emit a k-packet as TARDIS only works with photon packets.
Next we calculate the radidative
rates using equation 10 in :cite:`2003A&A...403..261L`.

.. math::
    {\cal R}_{\textrm{upper}\rightarrow\textrm{lower}} &=
    A_{\textrm{upper}\rightarrow\textrm{lower}}\beta_\textrm{Sobolev}n_\textrm{upper}\\
    {\cal R}_{\textrm{lower}\rightarrow\textrm{upper}} &=
    (B_{\textrm{lower}\rightarrow\textrm{upper}}n_\textrm{lower}-
    B_{\textrm{upper}\rightarrow\textrm{lower}}n_\textrm{upper})
    \beta_\textrm{Sobolev} J_{\nu}^{b}\\,

with :math:`\beta_\textrm{Sobolev} = \frac{1}{\tau_\textrm{Sobolev}}(1-e^{-\tau_\textrm{Sobolev}})` .

using the Einstein coefficients

.. math::
    A_{\textrm{upper}\rightarrow\textrm{lower}} &= \frac{8 \nu^2 \pi^2 e^2}{m_e c^3}~
        \frac{g_\textrm{lower}}{g_\textrm{upper}}~f_{\textrm{lower}\rightarrow\textrm{upper}}\\
    A_{\textrm{upper}\rightarrow\textrm{lower}} &= \underbrace{\frac{4 \pi^2 e^2}{m_e c}}_{C_\textrm{Einstein}}~ \frac{2\nu^2}{c^2}
            \frac{g_\textrm{lower}}{g_\textrm{upper}}~f_{\textrm{lower}\rightarrow\textrm{upper}}\\
    B_{\textrm{lower}\rightarrow\textrm{upper}} &= \frac{4\pi^2 e^2}{m_e h\nu c}\,f_{\textrm{lower}\rightarrow\textrm{upper}}\\
    B_{\textrm{lower}\rightarrow\textrm{upper}} &= \underbrace{\frac{4 \pi^2 e^2}{m_e c}}_{C_\textrm{Einstein}}\frac{1}{h\nu} f_{\textrm{lower}\rightarrow\textrm{upper}}\\

    B_{\textrm{upper}\rightarrow\textrm{lower}} &= \frac{4\pi^2 e^2}{m_e h\nu c}\,f_{\textrm{lower}\rightarrow\textrm{upper}}\\
    B_{\textrm{upper}\rightarrow\textrm{lower}} &= \underbrace{\frac{4 \pi^2 e^2}{m_e c}}_{C_\textrm{Einstein}}\frac{1}{h\nu}\frac{g_\textrm{lower}}{g_\textrm{upper}}f_{\textrm{lower}\rightarrow\textrm{upper}}\\

we get

.. math::
    {\cal R}_{\textrm{upper}\rightarrow\textrm{lower}} &=
        C_\textrm{Einstein} \frac{2\nu^2}{c^2} \frac{g_\textrm{lower}}{g_\textrm{upper}}~f_{\textrm{lower}\rightarrow\textrm{upper}}
        \beta_\textrm{Sobolev}n_\textrm{upper}\\

    {\cal R}_{\textrm{lower}\rightarrow\textrm{upper}} &=
            C_\textrm{Einstein}\frac{1}{h\nu} f_{\textrm{lower}\rightarrow\textrm{upper}}
            (n_\textrm{lower}-\frac{g_\textrm{lower}}{g_\textrm{upper}}n_\textrm{upper})
                        \beta_\textrm{Sobolev} J_{\nu}^{b}\\

This results in the transition probabilities:

.. math::
    p_\textrm{emission down}&= C_\textrm{Einstein} \frac{2\nu^2}{c^2} \frac{g_\textrm{lower}}{g_\textrm{i}}~f_{\textrm{lower}\rightarrow\textrm{i}}
                                       \beta_\textrm{Sobolev}n_\textrm{i}\,(\epsilon_\textrm{i} - \epsilon_\textrm{lower}) / {\cal D}_{i}\\
    p_\textrm{internal down}&= C_\textrm{Einstein} \frac{2\nu^2}{c^2} \frac{g_\textrm{lower}}{g_\textrm{i}}~f_{\textrm{lower}\rightarrow\textrm{i}}
                                       \beta_\textrm{Sobolev}n_\textrm{i}\,\epsilon_\textrm{lower}/{\cal D}_{i}\\
    p_\textrm{internal up}&=C_\textrm{Einstein}\frac{1}{h\nu} f_{\textrm{i}\rightarrow\textrm{upper}}
                                        (n_\textrm{i}-\frac{g_\textrm{i}}{g_\textrm{upper}}n_\textrm{upper})
                                                    \beta_\textrm{Sobolev} J_{\nu}^{b}\,\epsilon_{i}/{\cal D}_{i}\\,

and as we will normalise the transition probabilities numerically later,  we can get rid of :math:`C_\textrm{Einstein}`,
:math:`\frac{1}{{\cal D}_i}` and number densities.

.. math::
    p_\textrm{emission down}&= \frac{2\nu^2}{c^2} \frac{g_\textrm{lower}}{g_\textrm{i}}~f_{\textrm{lower}\rightarrow\textrm{i}}
                                       \beta_\textrm{Sobolev}\,(\epsilon_\textrm{i} - \epsilon_\textrm{lower})\\
    p_\textrm{internal down}&=  \frac{2\nu^2}{c^2} \frac{g_\textrm{lower}}{g_\textrm{i}}~f_{\textrm{lower}\rightarrow\textrm{i}}
                                       \beta_\textrm{Sobolev}\,\epsilon_\textrm{lower}\\
    p_\textrm{internal up}&=\frac{1}{h\nu} f_{\textrm{i}\rightarrow\textrm{upper}}
                                        \underbrace{(1-\frac{g_\textrm{i}}{g_\textrm{upper}}\frac{n_\textrm{upper}}{n_i})}
                                        _\textrm{stimulated emission}
                                                    \beta_\textrm{Sobolev} J_{\nu}^{b}\,\epsilon_{i}\\,




There are two parts for each of the probabilities, one that is pre-computed by `~tardisatomic` and is in the HDF5 File,
and one that is computed during the plasma calculations:

.. math::
        p_\textrm{emission down}&= \underbrace{\frac{2\nu^2}{c^2} \frac{g_\textrm{lower}}{g_\textrm{i}}~f_{\textrm{lower}\rightarrow\textrm{i}}
                                           (\epsilon_\textrm{i} - \epsilon_\textrm{lower})}_\textrm{pre-computed}
                                           \,\beta_\textrm{Sobolev}\\
        p_\textrm{internal down} &= \underbrace{\frac{2\nu^2}{c^2} \frac{g_\textrm{lower}}{g_\textrm{i}}~f_{\textrm{lower}\rightarrow\textrm{i}}
                                           \epsilon_\textrm{lower}}_\textrm{pre-computed}\,\beta_\textrm{Sobolev}\\
        p_\textrm{internal up} &= \underbrace{\frac{1}{h\nu} f_{\textrm{i}\rightarrow\textrm{upper}}}_\textrm{pre-computed}
                                                        \beta_\textrm{Sobolev} J_{\nu}^{b}\,
                                                        (1-\frac{g_\textrm{i}}{g_\textrm{upper}}\frac{n_\textrm{upper}}{n_i})
                                                        \,\epsilon_{i}.


.. automodapi:: tardisatomic.macro_atom_transition


Collisional Excitation
----------------------

Let us consider an ionized plasma, where atomic collisions with charged particles are the dominant process due to the
long range of the Coulomb interactions. In such a plasma, the collision rate is proportional to the flux of charged
particle. The high velocity of the electrons in comparison to those of ions allows us to consider only electrons for the
flux.

The first process described here is the Collisional Excitation where the translational energy of an electron collision
is converted into the internal energy of an atom/ion. As in most collisional process, Collisional Excitation
is characterised by a cross section. Let :math:`\sigma_{i,j,k \rightarrow i',j,k} (v)` be the cross section for
Collisional Excitation from the level :math:`i \rightarrow i'` due to an electron moving at :math:`v`.
Using the Maxwellian velocity distribution :maht:`f(v)` the number of Collisional Excitation  from :math:`i \rightarrow 1'` is

.. maht::
    n_i C_{i \rightarrow i'} = n_i n_e \int_{v_o}^{\infty} \sigma_{i \rightarrow i'} (v) f(v) dv.


In LTE we can obtain the number of Collisional Deexcitation via detailed balance :math:`C_{i\rightarrow i'} =
\frac{n_{i'}}{n_{i}} C_{i'\rightarrow i}`
In TARDIS we use the van Regemorter approximation (Equation 5-75 of Mihalas 1978) to approximate the collisional
cross-sections allowing the calculation of the collisional bound-bound rates for radiatively permitted transitions from
the oscillator strength. Hence,

.. math::
    c_{i \rightarrow i'} = n_e c_0 T^{1/2} 14.5 \frac{I_H}{h \nu_{i \rightarrow i'}} f_{i \rightarrow i'}
    \frac{h \nu_{i \rightarrow i'}}{k_b T}  exp\left(\frac{ - h \nu_{i \rightarrow i'}}{ k_b T} \right)
    \Zeta \left(\frac{h \nu_{i \rightarrow i'}}{k_b T}\right)

here :math:`c_0 = 5.46510E−11`, :math:`I_H = 13.6 \mathrm(eV)` is the ionization potential of hydrogen and :math:`\Zeta`
 is defined as

.. maht::
    \Zeta(u) = \max(a, 0.276 e^u E_1(u)).

E_1 is the Exponential integral and :math:`a` is defined as

.. maht::
    \[a =
        \left\{
            \begin{array}{lr}
                0.2 for n,l to n,l'\\
                0.7 for n,l to n'l'
            \end{array}
        \right.
    \]

Due to the lack of information about the main and spin quantum number for the transitions in the atomic dataset  TARDIS
assumes :math:`a = 0.5`.

Collisional Ionization
----------------------

Similary, we treat the collisional Ionization rate from the :maht:` i, j , k \rightarrow j+1 ,k`(in the following
:maht:` i \rightarrow j+1`  ). Hence the collisional Ionization rate is

.. math::
     n_i C_{i \rightarrow j+1} =  n_i n_e \int_{v_o}^{\infty} \sigma_{i \rightarrow j+1} (v) f(v) dv.

Here :maht:`\sigma_{i \rightarrow j+1}` is the photoionization cross-section and :maht:`f(v)` the Maxwellian velocity
distribution. In TARDIS we use the Seaton approximation (Equation 5-79 of Mihalas 1978) to obtain the collisional
bound-free rates

.. math::
    c_{i \rightarrow j+1} = n_e 1.55E13 T^{-0.5} b \frac{\sigma_{i \rightarrow j+1}(\nu_{i})}{h \nu_i /k_b /T} \exp\left(\frac{ - h \nu_{i}}{ k_b T} \right),

here :math:`b` is 0.1, 0.2 or 0.3 if the charge on the ionization state of level :math:`i` is 0, 1 or :maht:`\ge 2`.

Cooling via Free Free
---------------------

Cooling via Free Free is the conversion of heat into radiation energy by free free emission. Following Osterbrock 1974
the rate is given as

.. math::
    C = 1.426E-27 T^{-0.5} n_{0,j,k} n_e

Cooling via Free-Bound
----------------------
Cooling via Free-Bound describes the conversion of heat into radiant energy. The rate at which spontaneous recombination
transfers ionization and thermal energy into radiant energy is
:maht:`\alpha^{E,sp}_{i,j,k} h \nu_{i,j,k \rightarrow j+1,k} n_{j+1} n_e`
where
.. maht::
    \alpha^{E,sp}_i = 4 \pi  \int_{\nu_i}^{\infty} \frac{a_{i\kappa}(\nu)}{h \nu_i} \frac{2 h \nu^3}{c^2} e^{\frac{-h\nu}{k T}} d \nu,


.. matht::
    \alpha^{E,sp}_i = \frac{ 8 a_i k_B T  \nu_i^2 }{c^2 h} e^{frac{\nu_i h}{k_B T}}.


However, during such a recombination the ionization energy :maht:`h \nu_{i,j,k \rightarrow j+1,k}` is released as
radiant energy. The corresponding rate is given as

.. math::
        \alpha^{sp}_i = 4 \pi  \int_{\nu}^{\infty} \frac{a_{ik}(\nu)}{h\nu} \frac{2 h \nu^3}{c^2} e^{-h \nu / kt},


Subtracting this from the :maht:`\alpha^{E,sp}_i` and summing over all bound states, we get the rate for Cooling via Free-Bound as

.. math::
    C^{fb} = \nu_{i,j,k \rightarrow j+1,k} n_{j+1} n_e \sum_{i=0}^{I} \left(\alpha^{E,sp}_i - \alpha^{sp}_i \right) h \nu_{i,j,k \rightarrow j+1,k}.


Implementation
--------------

To speed up the plasma calculations during the runtime of TARDIS most of the equations mention above are spited up
into pre-computed coefficients which depend only on the atomic data and temperature. The coefficients are computed in tardisatomic and stored in the atomc HD5 file.




Collisional Excitation
----------------------

The Collisional Excitation in the van Regemorter approximation are pre computed in tardisatomic on a temperature grid.


.. math::
    c_{i \rightarrow i'} = n_e \underbrace{ c_0 T^{1/2} 14.5 \frac{I_H}{h \nu_{i \rightarrow i'}} f_{i \rightarrow i'}
    \frac{h \nu_{i \rightarrow i'}}{k_b T}  exp\left(\frac{ - h \nu_{i \rightarrow i'}}{ k_b T} \right)
    \Zeta \left(\frac{h \nu_{i \rightarrow i'}}{k_b T}\right)}_{pre-computed}

The pre-computed Collisional Excitation coefficients are stored in XXXXXXXX

    +--------------------+-----------------------+-----------------------+-------------------------------+-------------------------------+-------------------+---------------------------------------------------+
    |    atomic_number   |   source_ion_number   |   source_level_number |   destination_ion_number      |   destination_level_number    |   C_ul_conversion |   pre-computed coefficients from T[0] to T[-1]    |
    +--------------------+-----------------------+-----------------------+-------------------------------+-------------------------------+-------------------+---------------------------------------------------+
    |   int              |  int                  |  int                  |      int                      |   int                         |  int              |  double                                           |
    +--------------------+-----------------------+-----------------------+-------------------------------+-------------------------------+-------------------+---------------------------------------------------+

In TARDIS the pre-computed Collisional Excitation coefficients are then multiplied with the plasma depended factors, here


.. math:: c_{i \rightarrow i'} = n_e (\mahtrm{pre-computed coefficient})



Collisional Ionization
----------------------


The Collisional Ionization in the Seaton approximation are pre computed in tardisatomic on a temperature grid.

.. math::
    c_{i \rightarrow j+1} = n_e \underbrace{ 1.55E13 T^{-0.5} b \frac{\sigma_{i \rightarrow j+1}(\nu_{i})}{h \nu_i /k_b /T} \exp\left(\frac{ - h \nu_{i}}{ k_b T} \right)}_{pre-computed}


The pre-computed Collisional Excitation coefficients are stored in XXXXXXXX

    +--------------------+-----------------------+-----------------------+-------------------------------+-------------------------------+-------------------+---------------------------------------------------+
    |    atomic_number   |   source_ion_number   |   source_level_number |   destination_ion_number      |   destination_level_number    |   C_ul_conversion |   pre-computed coefficients from T[0] to T[-1]    |
    +--------------------+-----------------------+-----------------------+-------------------------------+-------------------------------+-------------------+---------------------------------------------------+
    |   int              |  int                  |  int                  |      int                      |   int                         |  int              |  double                                           |
    +--------------------+-----------------------+-----------------------+-------------------------------+-------------------------------+-------------------+---------------------------------------------------+

In TARDIS the pre-computed Collisional Excitation coefficients are then multiplied with the plasma depended factors, here

    .. math::
        c_{i \rightarrow j+1} = n_e (\mahtrm{pre-computed coefficient}).





Cooling via Free-Bound
----------------------

The computation of Cooling via Free-Bound is given as

.. math::
    C^{fb} = \nu_{i,j,k \rightarrow j+1,k} n_{j+1} n_e \sum_{i=0}^{I} \left(\underbrace{\alpha^{E,sp}_i}_{pre-computed} - \underbrace{\alpha^{sp}_i}_{pre-computed} \right) h \nu_{i,j,k \rightarrow j+1,k}.

This Equation can be spit up into to pre-computed factors  :math:`\alpha^{E,sp}_i` and :math:`\alpha^{sp}_i`.

To compute the values of :math:`\alpha^{E,sp}_i`  we use the analytical solution of the integral

.. math::
    alpha^{E,sp}_i = \\frac{ 8 a_i k_B T  \\nu_i^2 }{c^2 h} e^{frac{\\nu_i h}{k_B T}}.

Since there is no trivial analytical solution for the integral in the equation of :math:`\alpha^{sp}_i`, we use a numerical (Simpson rule) approach for the computation in tardisatomic.
The pre-computed values for both coefficients are computed on a on a temperature grid and stored in XXXXXX

    +--------------------+-----------------------+-----------------------+-------------------------------+-------------------------------+-------------------+---------------------------------------------------+---------------------------------------------------+
    |    atomic_number   |   source_ion_number   |   source_level_number |   destination_ion_number      |   destination_level_number    |   C_ul_conversion |   pre-computed coefficients SP from T[0] to T[-1] | pre-computed coefficients E,SP from T[0] to T[-1] |
    +--------------------+-----------------------+-----------------------+-------------------------------+-------------------------------+-------------------+---------------------------------------------------+---------------------------------------------------+
    |   int              |  int                  |  int                  |      int                      |   int                         |  int              |  double                                           | double
    +--------------------+-----------------------+-----------------------+-------------------------------+-------------------------------+-------------------+---------------------------------------------------+---------------------------------------------------+


In TARDIS the pre-computed coefficients are then used in combination  with the plasma depended factors to obtain the Free-Bound cooling rate, here


.. math::
    C^{fb} = \nu_{i,j,k \rightarrow j+1,k} n_{j+1} n_e \sum_{i=0}^{I} \left(\underbrace{\alpha^{E,sp}_i}_{pre-computed} - \underbrace{\alpha^{sp}_i}_{pre-computed} \right) h \nu_{i,j,k \rightarrow j+1,k}.









































