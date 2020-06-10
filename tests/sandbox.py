import selfiesv1 as sfv1


print(sfv1.encoder("CC1C(CCCC)1"))
print(sfv1.decoder(sfv1.encoder("CC1C(CCCC)1"), N_restrict=False))
