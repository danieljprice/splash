import subprocess

RENDERS = [
    ["splash", "2i_00000", "-r", "6",  "-p", "render_testing",                    "-dev", "log_rho_render.png"  ],
    ["splash", "2i_00000", "-r", "10", "-p", "render_testing",                    "-dev", "log_u_render.png"    ],
    ["splash", "2i_00000", "-r", "6",  "-p", "render_testing", "-vec", "7",       "-dev", "log_rho_v_render.png"],
    ["splash", "2i_00000", "-r", "6",  "-p", "render_testing", "--sink=0",        "-dev", "log_rho_sink0_render.png"],
    ["splash", "2i_00000", "-r", "6",  "-p", "render_testing", "--sink=1",        "-dev", "log_rho_sink1_render.png"],
    ["splash", "2i_00000", "-r", "6",  "-p", "render_testing", "--sink=2",        "-dev", "log_rho_sink2_render.png"],
]

for command in RENDERS:
    subprocess.run(command, check=True)