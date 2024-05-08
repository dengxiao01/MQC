import pyotp

key = 'TREQPZ4T5X4V7RFDIH7CXGGGVXKUWGR4'
totp = pyotp.TOTP(key)
print(totp.now())