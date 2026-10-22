# Draft: asking NCAR about the 3-D PBL scheme in ERF (22 Oct 2026)

Mike's suggestion (22 Oct 1:1): the 3-D PBL authors may already be putting it into ERF; ask them before building our own port.
Facts behind the email: ERF's MYNN documentation cites Juliano et al. 2022 and thanks Branko Kosović for discussions;
Kosović is a co-author of the ERF papers; the ERF source (v26.09 + development, checked 22 Oct) has MYNN2.5, MYNN-EDMF,
MRF and YSU but no 3-D PBL. Suggested recipients: Branko Kosović (NCAR), Timothy Juliano (NCAR); cc Mike and, if Mike
agrees, Eliot Quon / Ann Almgren (ERF leads). Send it from your MIT address after Mike has seen it.

---

**Subject:** 3-D PBL scheme in ERF — is a port under way?

Dear Dr. Kosović and Dr. Juliano,

I am a first-year PhD student in Michael Howland's group at MIT (cc'd). We run ERF alongside WRF for idealized
boundary-layer simulations on TACC's Stampede3, and my thesis direction is PBL-scheme development for the
mesoscale-to-LES range over complex terrain.

As a first step we would like to evaluate your three-dimensional PBL scheme (Juliano et al. 2022, MWR) against filtered
LES in ERF, using the same idealized cases we have already verified between WRF and ERF (GABLS1 stable, Shin & Dudhia
convective). ERF's documentation of its MYNN implementation cites your scheme, so before we start porting it ourselves
I wanted to ask:

1. Is a port of the 3-D PBL scheme into ERF already under way or planned on your side (or at the ERF team's)?
2. If not, would you be open to us doing it on a branch of the Howland Lab ERF fork, with your advice on the parts that
   differ from the WRF implementation (horizontal flux divergences on ERF's cell-centred grid, the surface-layer
   treatment of Eghdami et al. 2022)?
3. Is there a reference test case you would recommend for checking a port (we would start with the Perdigão
   configuration used in your WRF evaluations)?

Any pointers would be much appreciated. Thank you for your time.

Best regards,
Hyunsoo Na
PhD student, Civil and Environmental Engineering, MIT
Howland Lab

---

Before sending: (a) show Mike; (b) confirm the exact ERF version checked (26.09-538-ge50b1611bc98 on 22 Oct 2026);
(c) if the ERF team is cc'd, mention the `hsna` branch on github.com/Howland-Lab/ERF where our setup is recorded.
