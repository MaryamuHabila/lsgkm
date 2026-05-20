import re

code = open("run_tf_fimo_final.py").read()

code = re.sub(r'REF_FA.*', 'REF_FA="snp_ref_clean.fa"', code)
code = re.sub(r'ALT_FA.*', 'ALT_FA="snp_alt_clean.fa"', code)

open("run_tf_fimo_final.py","w").write(code)

print("FIMO paths fixed")
