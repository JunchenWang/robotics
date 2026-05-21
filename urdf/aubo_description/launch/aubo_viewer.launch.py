from launch import LaunchDescription
from launch.actions import DeclareLaunchArgument, OpaqueFunction
from launch.conditions import IfCondition, UnlessCondition
from launch.substitutions import Command, FindExecutable, LaunchConfiguration, PathJoinSubstitution
from launch_ros.actions import Node
from launch_ros.substitutions import FindPackageShare


def launch_setup(context, *args, **kwargs):
    robot_model = LaunchConfiguration("robot_model")
    transmission_hw_interface = LaunchConfiguration("transmission_hw_interface")
    use_joint_state_gui = LaunchConfiguration("use_joint_state_gui")
    rviz_config = LaunchConfiguration("rviz_config")

    robot_description_content = Command(
        [
            PathJoinSubstitution([FindExecutable(name="xacro")]),
            " ",
            PathJoinSubstitution(
                [FindPackageShare("aubo_description"), "urdf", "xacro", "inc", "aubo.xacro"]
            ),
            " ",
            "robot_model:=",
            robot_model,
            " ",
            "transmission_hw_interface:=",
            transmission_hw_interface,
        ]
    )
    robot_description = {"robot_description": robot_description_content}

    static_tf_node = Node(
        package="tf2_ros",
        executable="static_transform_publisher",
        name="static_transform_publisher",
        output="log",
        arguments=["0", "0", "0", "0", "0", "0", "world", "base_link"],
    )

    robot_state_publisher_node = Node(
        package="robot_state_publisher",
        executable="robot_state_publisher",
        name="robot_state_publisher",
        output="both",
        parameters=[robot_description],
    )

    joint_state_publisher_node = Node(
        package="joint_state_publisher",
        executable="joint_state_publisher",
        name="joint_state_publisher",
        output="log",
        condition=UnlessCondition(use_joint_state_gui),
        parameters=[robot_description],
    )

    joint_state_publisher_gui_node = Node(
        package="joint_state_publisher_gui",
        executable="joint_state_publisher_gui",
        name="joint_state_publisher_gui",
        output="log",
        condition=IfCondition(use_joint_state_gui),
        parameters=[robot_description],
    )

    rviz_node = Node(
        package="rviz2",
        executable="rviz2",
        name="rviz2",
        output="log",
        arguments=["-d", rviz_config],
        parameters=[robot_description],
    )

    return [
        static_tf_node,
        robot_state_publisher_node,
        joint_state_publisher_node,
        joint_state_publisher_gui_node,
        rviz_node,
    ]


def generate_launch_description():
    declared_arguments = [
        DeclareLaunchArgument(
            "robot_model",
            default_value="aubo_i5",
            description="Robot model to visualize, for example aubo_i5 or aubo_C3.",
        ),
        DeclareLaunchArgument(
            "transmission_hw_interface",
            default_value="hardware_interface/PositionJointInterface",
            description="Hardware interface inserted into the generated robot description.",
        ),
        DeclareLaunchArgument(
            "use_joint_state_gui",
            default_value="true",
            description="Whether to start joint_state_publisher_gui instead of joint_state_publisher.",
        ),
        DeclareLaunchArgument(
            "rviz_config",
            default_value=PathJoinSubstitution(
                [FindPackageShare("aubo_description"), "rviz", "view_robot.rviz"]
            ),
            description="Absolute path to the RViz config file.",
        ),
    ]

    return LaunchDescription(declared_arguments + [OpaqueFunction(function=launch_setup)])
