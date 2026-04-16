
energies = [f"lgE_{16.0 + i/10:.1f}" for i in range(21)] # 16.0 to 18.0
zeniths = [f"{i/10:.1f}" for i in range(10)]           # 0.0 to 0.9
primaries = ["proton", "iron", "helium", "oxygen"]               # Add others as needed


# with open("params.txt", "w") as f:
#     for p in primaries:
#         for e in energies:
#             for z in zeniths:
#                 f.write(f"{p} {e} {z}\n")



with open("params_plot.txt", "w") as f:
    for p in primaries:
        for e in energies:
            f.write(f"{p} {e}\n")