import matplotlib.pyplot as plt

from synaptus import (
    MultilayerCartesianPulseEchoData,
    MultilayerOmegaKMigration,
    PhaseShiftMigration,
    load_mat_dataset,
    plot_us_image,
)


def test_dataset_class_2d():
    mat_dataset = load_mat_dataset("LineScan2D_WireTargets.mat")

    dataset = MultilayerCartesianPulseEchoData(
        raw_data=mat_dataset.raw_data,
        fs=mat_dataset.fs,
        f_low=0.4e6,
        f_high=2.5e6,
        x_step=mat_dataset.x_step,  # type:ignore
        t_delay=mat_dataset.t_delay,
        wave_velocities=mat_dataset.wave_vel,
    )

    print(vars(dataset))


def test_phase_shift_migration_2d():
    """Test the PhaseShiftMigration class with 2D data."""
    mat_dataset = load_mat_dataset("LineScan2D_WireTargets.mat")

    psm = PhaseShiftMigration(
        raw_data=mat_dataset.raw_data,
        fs=mat_dataset.fs,
        f_low=0.4e6,
        f_high=2.5e6,
        x_step=mat_dataset.x_step,  # type:ignore
        t_delay=mat_dataset.t_delay,
        wave_velocities=mat_dataset.wave_vel,
    )
    images, z_vecs = psm.phase_shift_migrate()
    for image, z_vec in zip(images, z_vecs):
        plot_us_image(image, y_val=z_vec, title="PSM focused image, 2D, single layer")


def test_mulok_2d():
    """Test the MultilayerOmegaKMigration class with 2D data."""

    mat_dataset = load_mat_dataset("LineScan2D_WireTargets.mat")

    mulok = MultilayerOmegaKMigration(
        raw_data=mat_dataset.raw_data,
        fs=mat_dataset.fs,
        f_low=0.4e6,
        f_high=2.5e6,
        x_step=mat_dataset.x_step,  # type:ignore
        t_delay=mat_dataset.t_delay,
        wave_velocities=mat_dataset.wave_vel,
    )
    images, z_vecs = mulok.mulok_migrate()
    for image, z_vec in zip(images, z_vecs):
        plot_us_image(image, y_val=z_vec, title="MULOK focused image, 2D, single layer")


def test_psm_3d():
    """Test the PhaseShiftMigration class with 3D data."""
    mat_dataset = load_mat_dataset("PlaneScan3D_PlexiAluFBH.mat")

    psm = PhaseShiftMigration(
        raw_data=mat_dataset.raw_data,
        fs=mat_dataset.fs,
        f_low=mat_dataset.f_low,
        f_high=mat_dataset.f_high,
        x_step=mat_dataset.x_step,  # type:ignore
        y_step=mat_dataset.y_step,  # type:ignore
        t_delay=mat_dataset.t_delay,
        wave_velocities=mat_dataset.wave_vel,
        layer_thicknesses=mat_dataset.layer_thick,
    )

    # Perform migration (focusing)
    images, z_vecs = psm.phase_shift_migrate()

    # Create C-scan images plexiglas and aluminium by depth gating
    z_ind_plexi = (z_vecs[1] >= 0.068) & (z_vecs[1] <= 0.092)  # Depth gating, plexiglas layer
    z_ind_alu = (z_vecs[2] >= 0.1) & (z_vecs[2] <= 0.14)  # Depth gating, aluminium layer

    cscan_im_plexi = images[1][z_ind_plexi, :, :].max(axis=0)
    cscan_im_alu = images[2][z_ind_alu, :, :].max(axis=0)

    # Plot C-scan images
    fig, ax = plt.subplots(1, 2, figsize=(5, 5))
    ax[0].imshow(cscan_im_plexi)
    ax[1].imshow(cscan_im_alu)
    ax[0].set_title("C-scan, PSM, plexiglas")
    ax[1].set_title("C-scan, PSM, aluminium")
    plt.show()


def test_mulok_3d():
    """Test the MultilayerOmegaKMigration class with 3D data."""
    mat_dataset = load_mat_dataset("PlaneScan3D_PlexiAluFBH.mat")

    mulok = MultilayerOmegaKMigration(
        raw_data=mat_dataset.raw_data,
        fs=mat_dataset.fs,
        f_low=mat_dataset.f_low,
        f_high=mat_dataset.f_high,
        x_step=mat_dataset.x_step,  # type:ignore
        y_step=mat_dataset.y_step,  # type:ignore
        t_delay=mat_dataset.t_delay,
        wave_velocities=mat_dataset.wave_vel,
        layer_thicknesses=mat_dataset.layer_thick,
    )

    # Perform migration (focusing)
    images, z_vecs = mulok.mulok_migrate()

    # Create C-scan images plexiglas and aluminium by depth gating
    z_ind_plexi = (z_vecs[1] >= 0.068) & (z_vecs[1] <= 0.092)  # Depth gating, plexiglas layer
    z_ind_alu = (z_vecs[2] >= 0.1) & (z_vecs[2] <= 0.14)  # Depth gating, aluminium layer

    cscan_im_plexi = images[1][z_ind_plexi, :, :].max(axis=0)
    cscan_im_alu = images[2][z_ind_alu, :, :].max(axis=0)

    # Plot C-scan images
    fig, ax = plt.subplots(1, 2, figsize=(5, 5))
    ax[0].imshow(cscan_im_plexi)
    ax[1].imshow(cscan_im_alu)
    ax[0].set_title("C-scan, MULOK, plexiglas")
    ax[1].set_title("C-scan, MULOK, aluminium")
    plt.show()
