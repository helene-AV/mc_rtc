#include <mc_rtc/GUIState.h>
#include <mc_control/ControllerServer.h>

struct DummyProvider
{
  double value = 42.0;
  Eigen::Vector3d point = Eigen::Vector3d(0., 1., 2.);
};

struct TestServer
{
  TestServer();

  void publish();

  mc_control::ControllerServer server {1.0, 1.0, {"ipc:///tmp/mc_rtc_pub.ipc"}, {"ipc:///tmp/mc_rtc_rep.ipc"}};
  DummyProvider provider;
  mc_rtc::gui::StateBuilder builder;
  bool toggle_ = true;
};

TestServer::TestServer()
{
  builder.addElement({"dummy", "provider"}, mc_rtc::gui::Label("value",
                                                               [this](){ return provider.value; }));
  builder.addElement({"dummy", "provider"}, mc_rtc::gui::ArrayLabel("point",
                                                                    [this](){ return provider.point; }));
  builder.addElement({"Button example"}, mc_rtc::gui::Button("Push me",
                                                             [](){ LOG_INFO("Button pushed") }));
  builder.addElement({"Toggle example"}, mc_rtc::gui::Toggle("Toggle",
                                                             [this](){ return toggle_;},
                                                             [this](){ toggle_ = !toggle_; }));
}

void TestServer::publish()
{
  server.handle_requests(builder);
  server.publish(builder);
}

int main()
{
  TestServer server;
  while(1)
  {
    server.publish();
    usleep(50000);
  }
  return 0;
}
