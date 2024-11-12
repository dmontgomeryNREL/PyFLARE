Group Contribution Method
=========================

The **Group Contribution Method (GCM)**, as developed by Constantinou and Gani\ :footcite:p:`constantinou_new_1994,constantinou_estimation_1995` in the mid-1990s , provides a systematic approach for estimating the thermodynamic properties of pure organic compounds. This method decomposes molecules into structural groups, each contributing to a target property based on predefined group values. By summing these contributions, the GCM accurately predicts essential properties, including the acentric factor, normal boiling point, liquid molar volume at standard conditions (298 K) and more. This predictive capability is particularly useful for complex mixtures or novel compounds, where experimental thermodynamic data may be unavailable.  One example is sustainable aviation fuels (SAFs), which often have limited thermodynamic data available. GCM provides SAF developers with a means to estimate these critical properties without extensive physical testing, thereby aiding in the identification of promising fuel compositions before committing to large-scale production.

This GCM-Python code was developed with SAF research in mind. It builds on `Pavan B. Govindaraju's Matlab implementation <https://github.com/gpavanb-old/GroupContribution>`_ and has been expanded to include additional thermodynamic properties and mixture properties, calculated from the GCM values of individual compounds within the mixture. Table INSERT TABLE REFERENCE HERE outlines the properties for the *i-th* compound in a mixture, which depend on the *k-th* first-order and *j-th* second-order group contributions.

Table of GCM properties
-----------------------

.. table:: Props
   :widths: 20 10 30 20 40

   ==========================  =====================  ===========================================  ====================  ===========================================================
   Property                    Units                  Group Contributions                          Units                 Description
   ==========================  =====================  ===========================================  ====================  ===========================================================
   :math:`M_{w,i}`             kg/mol                 :math:`m_{w1k}`, :math:`m_{w2j}`             g/mol                 Molecular weight.
   :math:`T_{c,i}`             K                      :math:`t_{c1k}`, :math:`t_{c2j}`             1                     Critical temperature\ :footcite:p:`constantinou_new_1994`.
   :math:`p_{c,i}`             Pa                     :math:`p_{c1k}`, :math:`p_{c2j}`             bar\ :sup:`-0.5`      Critical pressure\ :footcite:p:`constantinou_new_1994`.
   :math:`V_{c,i}`             m\ :sup:`3`\ /mol      :math:`v_{c1k}`, :math:`v_{c2j}`             m\ :sup:`3`\ /kmol     
   :math:`T_{b,i}`             K                      :math:`t_{b1k}`, :math:`t_{b2j}`             1              
   :math:`T_{m,i}`             K                      :math:`t_{m1k}`, :math:`t_{m2j}`             1              
   :math:`\Delta H_{f,i}`      J/mol                  :math:`h_{f1k}`, :math:`h_{f2j}`             kJ/mol         
   :math:`\Delta G_{f,i}`      J/mol                  :math:`g_{f1k}`, :math:`g_{f2j}`             kJ/mol         
   :math:`\Delta H_{v,stp,i}`  J/mol                  :math:`h_{v1k}`, :math:`h_{v2j}`             kJ/mol         
   :math:`\omega_i`            1                      :math:`\omega_{1k}`, :math:`\omega_{2j}`     1              
   :math:`V_{m,stp,i}`         m\ :sup:`3`\ /mol      :math:`v_{m1k}`, :math:`v_{m2j}`             m\ :sup:`3`\ /kmol     
   :math:`C_{p,stp,i}`         J/mol/K                :math:`C_{pA1_k}`, :math:`C_{pA2_k}`, \      J/mol/K     
                                                      :math:`C_{pB1_k}`, :math:`C_{pB2_k}`, \
                                                      :math:`C_{pC1_k}`, :math:`C_{pC2_k}`  
   ==========================  =====================  ===========================================  ====================  ===========================================================


GCM individual properties
-------------------------

Properties for the *i-th* compound in a mixture are defined as the sum of contributions from first- and second-order groups. The equations are presented below:

.. math::

   \begin{align*}
   M_{w,i} &= \left[\sum_{k = 1}^{N_{g_1}}\mathbf{N}_{ik}m_{w1k} + W \sum_{j = 1}^{N_{g_2}} \mathbf{M}_{ij} m_{w2k} \right] \times 10^{-3} \\
   T_{c,i} &= 181.28 \ln  \left[ \sum_{k=1}^{N_{g_1}} \mathbf{N}_{ik} t_{c1k} + W \sum_{j=1}^{N_{g_2}} \mathbf{M}_{ij} t_{c2j} \right]
   \end{align*}

Individual correlations
-----------------------

Mixture properties from GCM
---------------------------

.. footbibliography::