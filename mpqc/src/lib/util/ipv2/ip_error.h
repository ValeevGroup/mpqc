
#define IPE_OK            0  /* No problem. */
#define IPE_KEY_NOT_FOUND 1  /* The keyword was not found. */
#define IPE_OUT_OF_BOUNDS 2  /* An array subscript was out of bounds. */
#define IPE_MALLOC        3  /* Memory allocation failed. */
#define IPE_NOT_AN_ARRAY  4  /* Gave index for data which isn't an array */
#define IPE_NOT_A_SCALAR  5  /* Didn't give index for data which is an array */
#define IPE_TYPE          6  /* The datum is not of the appropiate type. */
#define IPE_HAS_NO_VALUE  7  /* The keyword has no value. */
#define IPE_VAL_NOT_EXPD  8  /* A value was not expected for the keyword. */

