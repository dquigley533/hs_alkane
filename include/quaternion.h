/* Header file for C-compatible variables/functions in quaternion.F90 */

/* Find a minimum arc quaternion that rotates v1 onto v2 */
void quat_get_minimum_arc(double unit_vector1[3], double unit_vector2[3], double quat[4]);

/* Find quaternion that rotates around axis by angle */
void quat_axis_angle_to_quat(double rotation_axis[3], double angle, double quat[4]);

/* Find the product of two quaternions */
void quat_product(double quaternion1[4], double quaternion2[4], double quat_product[4], int normalise);

/* Find the inverse quaternion */
void quat_inverse(double quaternion[4], double quat_inverse[4]);

/* Conjugate a quaternion with a vector */
void quat_conjugate_q_with_v(double quaternion[4], double vector[3], double vec_out[3]);
