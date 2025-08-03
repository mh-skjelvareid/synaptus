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
        plot_us_image(image, y_val=z_vec)


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
        plot_us_image(image, y_val=z_vec)
